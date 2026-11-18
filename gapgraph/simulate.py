"""Synthetic finishing problem with known truth, for demos and tests.

True genome:   A  R  B  R  C  D  E          (R = 3 kb repeat, two copies)
Assembly graph:
  * R collapsed -> tangle  A->R, B->R, R->B, R->C   (resolve with long reads)
  * C | D  gap with no graph link                  (soft-clips + long reads + Hi-C)
  * D | E  no link and E stored reverse-complemented (orientation from all sources)
  * contig A carries a chimeric 2 kb insert from C  (clip hotspots -> split)
Reads: Illumina 150 bp, PacBio HiFi, ONT; Hi-C pairs.  Alignments are written
directly from the known truth (no aligner needed), with soft-clips exactly as
minimap2 -Y / bwa would produce at contig ends.
"""
import os
import random
from .seqs import revcomp, write_fasta

PLATFORMS = {
    "illumina": dict(len=(150, 150), sub=0.002, ins=0.0002, dele=0.0002, depth=20),
    "hifi": dict(len=(10000, 20000), sub=0.002, ins=0.0005, dele=0.0005, depth=12),
    "ont": dict(len=(8000, 30000), sub=0.03, ins=0.015, dele=0.02, depth=15),
}


def _rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def build_truth(rng, scale=1.0):
    sz = lambda n: int(n * scale)
    parts = {k: _rand_seq(rng, sz(n)) for k, n in
             dict(A=25000, R=3000, B=30000, C=20000, D=25000, E=20000).items()}
    order = ["A", "R", "B", "R", "C", "D", "E"]
    genome, spans, pos = [], [], 0
    for k in order:
        genome.append(parts[k])
        spans.append((k, pos, pos + len(parts[k])))
        pos += len(parts[k])
    genome = "".join(genome)
    span = {}
    for k, s, e in spans:
        span.setdefault(k, []).append((s, e))
    # contig blocks: (contig, contig_start, genome_start, length, orient)
    a0, a1 = span["A"][0]
    c0 = span["C"][0][0]
    ins = (c0 + sz(5000), sz(2000))
    cut = sz(12000)
    contigs = {
        "A": genome[a0:a0 + cut] + genome[ins[0]:ins[0] + ins[1]] + genome[a0 + cut:a1],
        "R": parts["R"], "B": parts["B"], "C": parts["C"], "D": parts["D"],
        "E": revcomp(parts["E"]),
    }
    blocks = [("A", 0, a0, cut, "+"), ("A", cut + ins[1], a0 + cut, a1 - a0 - cut, "+")]
    for s, e in span["R"]:
        blocks.append(("R", 0, s, e - s, "+"))
    for k in "BCD":
        s, e = span[k][0]
        blocks.append((k, 0, s, e - s, "+"))
    s, e = span["E"][0]
    blocks.append(("E", 0, s, e - s, "-"))
    links = [("A", "+", "R", "+"), ("R", "+", "B", "+"), ("B", "+", "R", "+"), ("R", "+", "C", "+")]
    truth = {"breakpoints_in_A": [cut, cut + ins[1]],
             "joins": [["C:R", "D:L"], ["D:R", "E:R"]],
             "repeat_paths": [["A+", "R+", "B+"], ["B+", "R+", "C+"]]}
    return genome, contigs, blocks, links, truth


def _simulate_read(rng, genome, gs, ge, p):
    """ops: list of (op, gpos, base) with op in M,I,D."""
    ops = []
    for g in range(gs, ge):
        r = rng.random()
        if r < p["ins"]:
            ops.append(("I", g, rng.choice("ACGT")))
        r = rng.random()
        if r < p["dele"]:
            ops.append(("D", g, ""))
            continue
        b = genome[g]
        if rng.random() < p["sub"]:
            b = rng.choice([x for x in "ACGT" if x != b])
        ops.append(("M", g, b))
    return ops


def _block_alignment(ops, blk, qname, rev_read, mapq=60):
    contig, cstart, gs, glen, orient = blk
    ge = gs + glen
    qbases = [b for op, _, b in ops if op != "D"]
    qlen = len(qbases)
    # query index of each op
    qi, idx = 0, []
    for op, g, b in ops:
        idx.append(qi)
        if op != "D":
            qi += 1
    inside = [i for i, (op, g, b) in enumerate(ops) if gs <= g < ge]
    if not inside:
        return None
    lo, hi = inside[0], inside[-1]
    while lo <= hi and ops[lo][0] != "M":
        lo += 1
    while hi >= lo and ops[hi][0] != "M":
        hi -= 1
    if hi - lo < 30:
        return None
    seg = ops[lo:hi + 1]
    lclip = idx[lo]
    rclip = qlen - (idx[hi] + 1)
    seq = "".join(qbases)
    if orient == "+":
        pos = cstart + (seg[0][1] - gs)
        cig = _cigar(seg)
    else:
        pos = cstart + (ge - 1 - seg[-1][1])
        cig = _cigar(seg[::-1])
        seq = revcomp(seq)
        lclip, rclip = rclip, lclip
    cigar = (f"{lclip}S" if lclip else "") + cig + (f"{rclip}S" if rclip else "")
    reverse = rev_read != (orient == "-")
    return {"qname": qname, "ref": contig, "pos": pos, "cigar": cigar, "seq": seq,
            "reverse": reverse, "aligned": hi - lo, "mapq": mapq,
            "qs": lclip if not reverse else rclip, "qe": qlen - (rclip if not reverse else lclip),
            "qlen": qlen, "rend": pos + sum(1 for op, _, _ in seg if op != "I")}


def _cigar(seg):
    out, last, n = [], None, 0
    for op, _, _ in seg:
        if op == last:
            n += 1
        else:
            if last:
                out.append(f"{n}{last}")
            last, n = op, 1
    out.append(f"{n}{last}")
    return "".join(out)


def simulate(outdir, seed=7, scale=1.0, hic_pairs=40000):
    rng = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)
    genome, contigs, blocks, links, truth = build_truth(rng, scale)
    G = len(genome)
    write_fasta([("genome", genome)], os.path.join(outdir, "truth.fasta"))
    write_fasta(list(contigs.items()), os.path.join(outdir, "contigs.fasta"))
    with open(os.path.join(outdir, "assembly.gfa"), "w") as fh:
        fh.write("H\tVN:Z:1.0\n")
        for n, s in contigs.items():
            fh.write(f"S\t{n}\t{s}\tLN:i:{len(s)}\n")
        for a, ao, b, bo in links:
            fh.write(f"L\t{a}\t{ao}\t{b}\t{bo}\t0M\n")
    files = {}
    for plat, p in PLATFORMS.items():
        n_reads = int(G * p["depth"] / sum(p["len"]) * 2)
        sam_path = os.path.join(outdir, f"{plat}.sam")
        with open(sam_path, "w") as fh:
            for n, s in contigs.items():
                fh.write(f"@SQ\tSN:{n}\tLN:{len(s)}\n")
            for r in range(n_reads):
                ln = rng.randint(*p["len"])
                gs = rng.randint(0, max(0, G - ln))
                ops = _simulate_read(rng, genome, gs, gs + ln, p)
                rev = rng.random() < 0.5
                qname = f"{plat}_{r}"
                alns = [a for a in (_block_alignment(ops, b, qname, rev) for b in blocks) if a]
                alns.sort(key=lambda a: -a["aligned"])
                for i, a in enumerate(alns):
                    flag = (16 if a["reverse"] else 0) | (2048 if i else 0)
                    fh.write(f"{qname}\t{flag}\t{a['ref']}\t{a['pos'] + 1}\t{a['mapq']}\t{a['cigar']}\t*\t0\t0\t{a['seq']}\t*\n")
        files[plat] = sam_path
    # Hi-C: power-law distance, 5% random noise
    pairs_path = os.path.join(outdir, "hic.pairs")
    gmap = []
    for blk in blocks:
        gmap.append(blk)

    def to_contig(g):
        for c, cs, gs, gl, o in gmap:
            if gs <= g < gs + gl:
                return c, (cs + g - gs) if o == "+" else (cs + gs + gl - 1 - g)
        return None

    with open(pairs_path, "w") as fh:
        fh.write("## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n")
        n = 0
        while n < hic_pairs:
            g1 = rng.randrange(G)
            if rng.random() < 0.05:
                g2 = rng.randrange(G)
            else:
                d = int(10 ** rng.uniform(3, 5.3))
                g2 = g1 + d if rng.random() < 0.5 else g1 - d
                if not 0 <= g2 < G:
                    continue
            m1, m2 = to_contig(g1), to_contig(g2)
            if not m1 or not m2:
                continue
            fh.write(f"hic{n}\t{m1[0]}\t{m1[1] + 1}\t{m2[0]}\t{m2[1] + 1}\t+\t-\n")
            n += 1
    import json
    with open(os.path.join(outdir, "truth.json"), "w") as fh:
        json.dump(truth, fh, indent=1)
    return {"gfa": os.path.join(outdir, "assembly.gfa"), "sam": files, "pairs": pairs_path, "truth": truth}


# ---------------------------------------------------------------------------
# circular bacterial genome
# ---------------------------------------------------------------------------
CIRC_PLATFORMS = {
    "illumina": dict(len=(150, 150), sub=0.002, ins=0.0002, dele=0.0002, depth=30),
    "ont": dict(len=(3000, 12000), sub=0.03, ins=0.015, dele=0.02, depth=30),
}


def _mutate(rng, s, rate):
    return "".join(rng.choice([x for x in "ACGT" if x != c]) if rng.random() < rate else c for c in s)


def simulate_circular(outdir, seed=11, scale=0.5):
    """Chromosome (dnaA at base 1 in the truth) + plasmid.

    Assembly as an assembler might leave it:
      chr_1  forward, starting at 30% of the chromosome (not at dnaA)
      chr_2  reverse-complemented, overlaps chr_1 by 300 bp at the break and
             runs past the origin, duplicating 1 kb of chr_1's start
      plasmid  single contig with a 600 bp duplicated end overlap
    No GFA links: closure must come from reads.
    """
    rng = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)
    G = int(200000 * scale)
    P = int(16000 * scale)
    dnaA = "ATG" + _rand_seq(rng, 1398) + "TAA"
    chrom = dnaA + _rand_seq(rng, G - len(dnaA))
    plasmid = _rand_seq(rng, P)
    o, s1, brk_ovl, ori_ovl, p_ovl = int(0.3 * G), int(0.45 * G), 300, 1000, 600
    x_start = o + s1 - brk_ovl
    X = chrom[x_start:] + chrom[:o + ori_ovl]
    contigs = {"chr_1": chrom[o:o + s1], "chr_2": revcomp(X), "plasmid": plasmid + plasmid[:p_ovl]}
    # blocks are contiguous in "unrolled" molecule coordinates (copies of the
    # circle laid end to end), so reads crossing the origin align as one piece
    chrom_blocks = [("chr_1", 0, o, s1, "+"), ("chr_2", 0, x_start, len(X), "-")]
    plasmid_blocks = [("plasmid", 0, 0, P + p_ovl, "+")]

    write_fasta([("chromosome", chrom), ("plasmid", plasmid)], os.path.join(outdir, "truth.fasta"))
    write_fasta(list(contigs.items()), os.path.join(outdir, "contigs.fasta"))
    write_fasta([("dnaA_related_strain", _mutate(rng, dnaA, 0.04))], os.path.join(outdir, "dnaA_related.fasta"))
    with open(os.path.join(outdir, "assembly.gfa"), "w") as fh:
        fh.write("H\tVN:Z:1.0\n")
        for n, s in contigs.items():
            fh.write(f"S\t{n}\t{s}\tLN:i:{len(s)}\n")

    for plat, p in CIRC_PLATFORMS.items():
        with open(os.path.join(outdir, f"{plat}.sam"), "w") as fh:
            for n, s in contigs.items():
                fh.write(f"@SQ\tSN:{n}\tLN:{len(s)}\n")
            r = 0
            for mol, blocks in ((chrom, chrom_blocks), (plasmid, plasmid_blocks)):
                M = len(mol)
                ext = mol * (2 + p["len"][1] // M)          # reads wrap the origin
                copies = len(ext) // M
                shifted = [(c, cs, gs + k * M, gl, ori) for k in range(-1, copies) for c, cs, gs, gl, ori in blocks]
                for _ in range(int(M * p["depth"] / (sum(p["len"]) / 2))):
                    ln = rng.randint(*p["len"])
                    gs = rng.randrange(M)
                    ops = _simulate_read(rng, ext, gs, gs + ln, p)
                    rev = rng.random() < 0.5
                    qname = f"{plat}_{r}"
                    r += 1
                    alns = [a for a in (_block_alignment(ops, b, qname, rev) for b in shifted) if a]
                    alns.sort(key=lambda a: -a["aligned"])
                    for i, a in enumerate(alns):
                        flag = (16 if a["reverse"] else 0) | (2048 if i else 0)
                        fh.write(f"{qname}\t{flag}\t{a['ref']}\t{a['pos'] + 1}\t{a['mapq']}\t{a['cigar']}"
                                 f"\t*\t0\t0\t{a['seq']}\t*\n")
    truth = {"chromosome_length": G, "plasmid_length": P,
             "breaks": "chr_1:R ~ chr_2:R overlap 300; chr_2:L ~ chr_1:L overlap 1000 (origin)",
             "plasmid_end_overlap": p_ovl, "dnaA": "chromosome base 1, + strand"}
    import json
    with open(os.path.join(outdir, "truth.json"), "w") as fh:
        json.dump(truth, fh, indent=1)
    return {"gfa": os.path.join(outdir, "assembly.gfa"),
            "sam": {pl: os.path.join(outdir, f"{pl}.sam") for pl in CIRC_PLATFORMS},
            "start": os.path.join(outdir, "dnaA_related.fasta"), "truth": truth}
