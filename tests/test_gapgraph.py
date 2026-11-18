import os
import pytest
from gapgraph.simulate import simulate
from gapgraph.session import Session
from gapgraph.seqs import read_fasta, revcomp
from gapgraph.gfa import AssemblyGraph, join_key
from gapgraph.alignments import parse_cigar, read_sam
from gapgraph.project import Project, ProjectError


@pytest.fixture(scope="module")
def sim(tmp_path_factory):
    d = str(tmp_path_factory.mktemp("sim"))
    r = simulate(d, scale=0.3, hic_pairs=12000, seed=3)
    S = Session(gfa=r["gfa"], alignments=[(r["sam"][p], p) for p in ("illumina", "hifi", "ont")],
                pairs=r["pairs"])
    return d, r, S


def test_revcomp():
    assert revcomp("AACGTN") == "NACGTT"


def test_link_ends_are_symmetric(tmp_path):
    g = tmp_path / "g.gfa"
    g.write_text("S\ta\tACGT\nS\tb\tGGCC\nL\ta\t+\tb\t-\t2M\n")
    G = AssemblyGraph.from_gfa(str(g))
    assert G.links[0].ends == (("a", "R"), ("b", "R"))
    assert G.links[0].overlap == 2
    assert set(G.dead_ends()) == {("a", "L"), ("b", "L")}


def test_cigar():
    assert parse_cigar("5S10M2I3D4H") == [(5, "S"), (10, "M"), (2, "I"), (3, "D"), (4, "H")]


def test_sam_consistent_with_contigs(sim):
    d, r, S = sim
    seqs = read_fasta(os.path.join(d, "contigs.fasta"))
    m = t = 0
    for a in read_sam(r["sam"]["hifi"]):
        q, p = a.lclip_soft, a.pos
        for n, op in a.cigar:
            if op == "M":
                m += sum(a.seq[q + i] == seqs[a.ref][p + i] for i in range(n))
                t += n
                q, p = q + n, p + n
            elif op == "I":
                q += n
            elif op == "D":
                p += n
        assert p == a.end
    assert m / t > 0.99


def test_true_joins_found_and_ranked(sim):
    _, _, S = sim
    top = {c["id"] for c in S.candidates[:10]}
    assert "C:R~D:L" in top
    assert "D:R~E:R" in top          # E is stored reverse-complemented
    c = S.cand_by_id["C:R~D:L"]
    assert c["longread_reads"] > 0 and c["softclip_reads"] > 0 and c["hic_contacts"] > 0
    assert not c["graph_link"]


def test_repeat_resolution_paths(sim):
    _, _, S = sim
    assert S.graph.tangles() == ["R"]
    vias = {c["id"]: c["via"] for c in S.candidates if c["via"]}
    assert vias.get("A:R~B:L") == ["R+"]
    assert vias.get("B:R~C:L") == ["R+"]


def test_misjoin_hotspots(sim):
    d, r, S = sim
    bp = set(r["truth"]["breakpoints_in_A"])
    found = {h["pos"] for h in S.hotspots if h["contig"] == "A"}
    assert bp <= found


def test_full_reconstruction(sim, tmp_path):
    d, r, S = sim
    S.project = Project(S.seqs)
    b1, b2 = r["truth"]["breakpoints_in_A"]
    S.op({"op": "split", "piece": "A", "pos": b1})
    S.op({"op": "split", "piece": f"A:{b1}-{len(S.seqs['A'])}", "pos": b2 - b1})
    S.op({"op": "join", "a": [f"A:0-{b1}", "R"], "b": [f"A:{b2}-{len(S.seqs['A'])}", "L"], "gap": 0})
    chim = [s["name"] for s in S.project.scaffolds if s["parts"][0]["piece"] == f"A:{b1}-{b2}"][0]
    S.op({"op": "discard", "scaffold": chim})
    S.op({"op": "auto"})
    S.project.export(str(tmp_path / "out"))
    genome = read_fasta(os.path.join(d, "truth.fasta"))["genome"]
    out = read_fasta(str(tmp_path / "out.fasta"))
    assert any(s == genome or revcomp(s) == genome for s in out.values())


def test_undo_and_errors():
    P = Project({"x": "A" * 50, "y": "C" * 50})
    P.join(("x", "R"), ("y", "R"), gap=10)
    assert len(P.scaffolds) == 1
    assert P.scaffold_seq(P.scaffolds[0]) == "A" * 50 + "N" * 10 + "G" * 50
    with pytest.raises(ProjectError):
        P.join(("x", "L"), ("y", "L"))   # same scaffold
    P.undo()
    assert len(P.scaffolds) == 2


def test_overlap_trim():
    P = Project({"x": "AAAACCGG", "y": "CCGGTTTT"})
    P.join(("x", "R"), ("y", "L"), gap=-4)
    assert P.scaffold_seq(P.scaffolds[0]) == "AAAACCGGTTTT"
