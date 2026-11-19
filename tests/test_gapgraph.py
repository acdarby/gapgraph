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
    d, r, _ = sim
    S = Session(gfa=r["gfa"], alignments=[(r["sam"][p], p) for p in ("illumina", "hifi", "ont")],
                pairs=r["pairs"])
    b1, b2 = r["truth"]["breakpoints_in_A"]
    LA = len(S.orig_seqs["A"])
    S.op({"op": "split", "contig": "A", "pos": b1})
    # the split is real: two new contigs, reads divided between them
    assert "A" not in S.seqs and {f"A:0-{b1}", f"A:{b1}-{LA}"} <= set(S.seqs)
    assert S.stores["illumina"].by_ref[f"A:{b1}-{LA}"]
    S.op({"op": "split", "contig": f"A:{b1}-{LA}", "pos": b2 - b1})
    S.op({"op": "join", "a": [f"A:0-{b1}", "R"], "b": [f"A:{b2}-{LA}", "L"], "gap": 0})
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
        P.join(("x", "R"), ("x", "R"))   # same end twice
    P.join(("x", "L"), ("y", "L"), gap=0)  # remaining ends -> circularise
    assert P.scaffolds[0]["circular"]
    P.undo()
    P.undo()
    assert len(P.scaffolds) == 2


def test_overlap_trim():
    P = Project({"x": "AAAACCGG", "y": "CCGGTTTT"})
    P.join(("x", "R"), ("y", "L"), gap=-4)
    assert P.scaffold_seq(P.scaffolds[0]) == "AAAACCGGTTTT"


# ------------------------------------------------------------ circular genomes
from gapgraph.simulate import simulate_circular
from gapgraph import circular


@pytest.fixture(scope="module")
def circ(tmp_path_factory):
    d = str(tmp_path_factory.mktemp("circ"))
    r = simulate_circular(d, scale=0.3, seed=5)
    return d, r


def _session(r, **kw):
    return Session(gfa=r["gfa"], alignments=[(r["sam"]["illumina"], "illumina"), (r["sam"]["ont"], "ont")], **kw)


def test_overlap_and_start_detection():
    import random
    rng = random.Random(2)
    g = "".join(rng.choice("ACGT") for _ in range(30000))
    assert circular.end_overlap(g + g[:777]) == 777
    assert circular.end_overlap(g) == 0
    gene = g[12000:13400]
    hit = circular.find_start(g, revcomp(gene))
    assert hit["strand"] == "-" and hit["pos"] == 13399
    # gene spanning the current origin
    rot = g[12500:] + g[:12500]
    assert circular.find_start(rot, gene)["pos"] == len(g) - 500


def test_wrapping_reads_are_evidence(circ):
    d, r = circ
    S = _session(r)
    c = S.cand_by_id["plasmid:L~plasmid:R"]
    assert c["longread_reads"] > 5 and c["softclip_reads"] > 5
    assert c["gap"] == -r["truth"]["plasmid_end_overlap"]   # snapped to exact overlap
    assert S.cand_by_id["chr_1:L~chr_2:L"]["gap"] == -1000
    assert S.cand_by_id["chr_1:R~chr_2:R"]["gap"] == -300


def test_bacterial_auto_finish_is_exact(circ, tmp_path):
    d, r = circ
    S = _session(r, start_fasta=r["start"])
    res = S.auto_scaffold()
    assert len(res["circularised"]) == 2 and res["rotated"]
    S.project.export(str(tmp_path / "fin"))
    truth = read_fasta(os.path.join(d, "truth.fasta"))
    out = read_fasta(str(tmp_path / "fin.fasta"))
    assert sorted(out.values(), key=len) == sorted(truth.values(), key=len)
    agp = (tmp_path / "fin.agp").read_text()
    assert "circular" in agp


def test_rotation_flip_consistency():
    seq = "ATGAAACCCGGGTTTAAACCCGGGTTTA" * 5
    P = Project({"c": seq + seq[:30]})
    P.circularise("scf_c", -30)
    P.set_origin("scf_c", 40, "+")
    s1 = P.scaffold_seq(P.scaffolds[0])
    P.flip("scf_c")
    assert P.scaffold_seq(P.scaffolds[0]) == s1   # same molecule, same start, same strand
    assert s1 == seq[40:] + seq[:40]
    P.flip("scf_c")
    circle = P.scaffold_seq(P.scaffolds[0], rotate=False)
    P.set_origin("scf_c", 50, "-")          # gene's first base is circle[50], on the minus strand
    out = P.scaffold_seq(P.scaffolds[0])
    rc = revcomp(circle)
    k = len(circle) - 1 - 50
    assert out == rc[k:] + rc[:k] and out[0] == revcomp(circle[50])


# ------------------------------------------- consensus, problems, pairs, join editor
from gapgraph.alignments import Alignment, parse_cigar, AlignmentStore
from gapgraph.consensus import Pileup
from gapgraph import pairs as pairs_mod, joinedit
from collections import Counter


def _aln(name, ref, pos, cigar, seq, platform="illumina", flag=0, mate_ref="", mate_pos=-1, tlen=0, mapq=60):
    cg = parse_cigar(cigar)
    rlen = sum(n for n, op in cg if op in "MD")
    return Alignment(name, ref, pos, pos + rlen, bool(flag & 16), mapq, cg, seq, 0, 0, 0, 0, False,
                     len(seq), platform, flag, mate_ref, mate_pos, tlen, "")


def test_consensus_finds_substitution_deletion_insertion():
    import random
    rng = random.Random(4)
    contig = "".join(rng.choice("ACGT") for _ in range(300))
    # truth differs from the contig: sub at 50, contig has an extra base at 150,
    # contig is missing a base before 250
    truth = contig[:50] + ("A" if contig[50] != "A" else "C") + contig[51:150] + contig[151:250] + "G" + contig[250:]
    alns = []
    for k in range(12):
        # reads = truth; express their alignment to the contig
        seq = truth
        cig = "150M1D99M1I50M"   # 0-149 match(50 is a mismatch), skip contig 150, ..., insert G before 250
        alns.append(_aln(f"r{k}", "c", 0, cig, seq))
    P = Pileup("c", len(contig)).add_alignments(alns)
    regs = P.regions(P.masks(contig, "contig", 3), merge=0)
    assert [r[0] for r in regs] == [50, 150, 250]
    assert P.qual[10] == 100 and P.consensus_string(150, 151) == "*"


def test_problem_search_finds_planted_error(circ):
    d, r = circ
    S = _session(r)
    s = S.seqs["chr_1"]
    S.seqs["chr_1"] = s[:2000] + ("A" if s[2000] != "A" else "C") + s[2001:]
    S._pileups.clear()
    assert S.problems("chr_1", "contig") == [[2000, 2001]]
    assert S.next_problem("chr_1", 0, "contig")["start"] == 2000
    assert S.next_problem("chr_1", 2500, "contig") is None
    assert S.next_problem("chr_1", 2500, "contig", direction=-1)["start"] == 2000


def test_pairs_insert_size_and_spanning(circ):
    d, r = circ
    S = _session(r)
    st = S.insert_stats["illumina"]
    assert 470 <= st["median"] <= 530
    c = S.cand_by_id["chr_1:R~chr_2:R"]
    assert c["pair_reads"] >= 10 and abs(c["pair_gap"] + 300) < 40
    t = S.templates("chr_1", 0, S.lengths["chr_1"])
    assert t["counts"].get("consistent", 0) > 100 and t["counts"].get("spanning", 0) > 0


def test_pair_classification():
    stats = {"illumina": {"lo": 300, "hi": 700, "median": 500, "sd": 50}}
    ok = _aln("p", "c", 100, "150M", "A" * 150, flag=1 | 32 | 64, mate_ref="c", mate_pos=450, tlen=500)
    assert pairs_mod.classify(ok, stats) == "consistent"
    same = _aln("p", "c", 100, "150M", "A" * 150, flag=1 | 64, mate_ref="c", mate_pos=450, tlen=500)
    assert pairs_mod.classify(same, stats) == "orientation"
    far = _aln("p", "c", 100, "150M", "A" * 150, flag=1 | 32 | 64, mate_ref="c", mate_pos=5000, tlen=5050)
    assert pairs_mod.classify(far, stats) == "size"
    other = _aln("p", "c", 100, "150M", "A" * 150, flag=1 | 32 | 64, mate_ref="d", mate_pos=10)
    assert pairs_mod.classify(other, stats) == "spanning"


def test_join_view_overlaps_agree_in_register(circ):
    d, r = circ
    S = _session(r)
    for cid in ("chr_1:R~chr_2:R", "chr_1:L~chr_2:L", "plasmid:L~plasmid:R"):
        g = S.cand_by_id[cid]["gap"]
        v = S.joinview({"id": cid}, g - 30, 30)
        cols = [(v["A"]["consensus"][x - v["start"]], v["B"]["consensus"][x - v["start"]]) for x in range(g, 0)]
        assert cols and all(a == b for a, b in cols), cid
        al = S.joinalign({"id": cid})
        assert al["identity"] == 100.0 and al["overlap"] == -g


def test_banded_align_indels():
    a = "ACGTACGTTAGCATCGATCGGATCGATTACG"
    b = a[:10] + a[11:20] + "T" + a[20:]
    x, y = joinedit.banded_align(a, b)
    assert x.replace("*", "") == a and y.replace("*", "") == b
    assert x.count("*") == 1 and y.count("*") == 1


# ------------------------------- working assembly: splits and joins are real edits
from gapgraph import realign
from gapgraph.working import cut_alignment, reverse_alignment


def test_reverse_and_cut_alignment():
    a = _aln("r", "c", 100, "5S20M2I10M3D15M4S", "ACGT" * 14)   # 56 query bases
    a = replace_soft(a, 5, 4)
    b = reverse_alignment(reverse_alignment(a, 1000), 1000)
    assert (b.pos, b.end, b.cigar, b.seq) == (a.pos, a.end, a.cigar, a.seq)
    r = reverse_alignment(a, 1000)
    assert r.pos == 1000 - a.end and r.lclip_soft == 4 and r.rclip_soft == 5
    c = cut_alignment(a, 0, 125)          # keep the first 25 ref bases
    assert c.end == 125 and c.pos == 100 and c.lclip_soft == 5
    assert c.rclip_soft == len(a.seq) - 5 - sum(n for n, op in c.cigar if op in "MI")


def replace_soft(a, l, r):
    from dataclasses import replace
    return replace(a, lclip=l, lclip_soft=l, rclip=r, rclip_soft=r)


def test_split_makes_real_contigs_and_undo_restores(circ):
    d, r = circ
    S = _session(r)
    res = S.op({"op": "split", "contig": "plasmid", "pos": 4000})
    assert res["gone"] == ["plasmid"] and len(res["new"]) == 2
    left = [w for w in res["new"] if w.startswith("plasmid:0-")][0]
    # reads that crossed the break are now soft-clipped at the new end
    L = S.lengths[left]
    clipped = [a for a in S.stores["ont"].fetch(left, L - 50, L) if a.rclip_soft > 50 and a.end >= L - 5]
    assert len(clipped) >= 3
    # ... and the evidence immediately proposes joining it back
    back = [c for c in S.candidates if {c["a"][0], c["b"][0]} == set(res["new"])]
    assert back and back[0]["longread_reads"] > 0
    S.op({"op": "undo"})
    assert "plasmid" in S.seqs


def test_join_makes_one_contig_with_reads_across(circ):
    d, r = circ
    S = _session(r)
    res = S.op({"op": "join", "id": "chr_1:R~chr_2:R"})
    w = res["new"][0]
    j = S.W.comps[w][1][0]
    reads = [a for st in S.stores.values() for a in st.fetch(w, j - 5, j + 5) if not a.supplementary]
    assert len({(a.qname, a.flag & 192) for a in reads}) == len(reads)        # no duplicates
    assert sum(1 for a in reads if a.pos < j - 50 and a.end > j + 50) >= 15   # reads span the join
    assert not [p for p in S.problems(w, "contig") if abs(p[0] - j) < 2000]
    st = Counter(pairs_mod.classify(a, S.insert_stats, S.other_contigs)
                 for a in S.stores["illumina"].fetch(w, j - 300, j + 300) if not a.supplementary)
    assert st.get("spanning", 0) == 0 and st["consistent"] > 20


def test_clips_align_to_partner_in_join_editor(circ):
    d, r = circ
    S = _session(r)
    v = S.joinview({"id": "chr_1:R~chr_2:R"}, -400, 300)
    s0, ref, m, t = v["start"], v["B"]["ref"], 0, 0
    for row in v["A"]["rows"]:
        for rd in row:
            for i, ch in enumerate(rd["s"]):
                x = rd["x"] + i
                if ch.islower() and ch != "*" and not any(u0 <= x < u1 for u0, u1 in rd.get("u", [])):
                    o = ref[x - s0] if 0 <= x - s0 < len(ref) else " "
                    if o != " ":
                        t += 1
                        m += ch.upper() == o
    assert t > 1000 and m / t > 0.95


def test_clip_consensus_from_noisy_reads():
    import random
    from gapgraph.simulate import _simulate_read, PLATFORMS
    rng = random.Random(9)
    truth = "".join(rng.choice("ACGT") for _ in range(1500))
    clips = []
    for _ in range(15):
        ops = _simulate_read(rng, truth, 0, rng.randint(900, 1500), PLATFORMS["ont"])
        clips.append("".join(b for op, g, b in ops if op != "D"))
    cons = realign.clip_consensus(clips)
    assert len(cons) > 800
    x, y = joinedit.banded_align(cons, truth[:len(cons)])
    assert sum(p == q for p, q in zip(x, y)) / len(x) > 0.99


# --------------------------------------------------------------------- foldbacks
from gapgraph import foldback


def test_end_foldback_geometry():
    import random
    rng = random.Random(2)
    g = "".join(rng.choice("ACGT") for _ in range(30000))
    hit = foldback.end_foldback(g + revcomp(g[-1500:]), "R")
    assert hit["fold"] == 30000 and hit["arm"] == 1500 and hit["identity"] == 1.0
    loop = "".join(rng.choice("ACGT") for _ in range(30))
    assert abs(foldback.end_foldback(g + loop + revcomp(g[-1500:]), "R")["fold"] - 30015) <= 1
    assert foldback.end_foldback(revcomp(g[:2000]) + g, "L")["fold"] == 2000
    assert foldback.end_foldback(g, "R") is None and foldback.end_foldback(g, "L") is None


def test_real_hairpin_with_spanning_reads_is_not_auto_removed():
    import random
    from gapgraph.alignments import AlignmentStore
    rng = random.Random(3)
    g = "".join(rng.choice("ACGT") for _ in range(8000))
    seq = g + revcomp(g[-1500:])
    fold = len(g)
    reads = [_aln(f"s{i}", "c", fold - 400 + i, "800M", seq[fold - 400 + i:fold + 400 + i]) for i in range(10)]
    res = foldback.contig_foldbacks({"c": seq}, {"ont": AlignmentStore().add(reads)})
    assert res and res[0]["spanning_reads"] == 10 and res[0]["recommend"] == "check"


@pytest.fixture(scope="module")
def fbsim(tmp_path_factory):
    d = str(tmp_path_factory.mktemp("fb"))
    r = simulate_circular(d, scale=0.4, seed=5, foldback_reads=0.08, foldback_contig=True)
    return d, r


def test_foldback_reads_detected_exactly_and_trimmed(fbsim):
    d, r = fbsim
    S = _session(r, start_fasta=r["start"])
    ev = S.fb_events["ont"]
    assert len(ev) == r["truth"]["foldback_reads"] > 0
    assert S.fb_events["illumina"] == []
    assert S.fb_stats["ont"]["at_contig_foldback"] > 0      # hairpin-arm reads not miscounted
    fb_names = {e["qname"] for e in ev}
    left = [a for a in S.orig_stores["ont"].all() if a.qname in fb_names]
    assert len(left) == len(fb_names)                        # folded arm alignment removed
    assert all(a.rclip_soft == 0 and a.rclip > 0 for a in left)  # folded arm hard-clipped
    S.set_foldback_mode("drop")
    assert not [a for a in S.orig_stores["ont"].all() if a.qname in fb_names]


def test_foldback_contig_end_removed_and_genome_exact(fbsim, tmp_path):
    d, r = fbsim
    S = _session(r, start_fasta=r["start"])
    t = r["truth"]["foldback_contig"]
    fb = S.contig_foldbacks
    assert len(fb) == 1 and (fb[0]["contig"], fb[0]["side"], fb[0]["fold"], fb[0]["arm"]) == \
        (t["contig"], t["side"], t["fold"], t["arm"])
    assert fb[0]["recommend"] == "remove"
    S.op({"op": "remove_foldback", "contig": "chr_1", "side": "R"})
    assert not S.contig_foldbacks and len(S.project.discarded) == 1
    S.op({"op": "undo"})                                     # one undo reverses split + discard
    assert "chr_1" in S.seqs and not S.project.discarded
    res = S.auto_scaffold()
    assert res["foldbacks_removed"] and len(res["circularised"]) == 2
    S.project.export(str(tmp_path / "fin"))
    truth = read_fasta(os.path.join(d, "truth.fasta"))
    out = read_fasta(str(tmp_path / "fin.fasta"))
    assert sorted(out.values(), key=len) == sorted(truth.values(), key=len)
