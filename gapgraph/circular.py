"""Circular replicons (bacterial chromosomes, plasmids, organelles).

Evidence that a contig (or scaffold) closes on itself:
  * GFA self-link  X+ -> X+  (Flye, hifiasm, Unicycler, Verkko write these)
  * assembler naming: hifiasm 'ptg000001c' / 'utg...c' suffix
  * reads that run off the R end and continue at the L end: soft-clips and
    long reads (primary + supplementary) wrapping the origin
  * a duplicated suffix==prefix overlap (Canu keeps one; others sometimes)
Finishing then = trim the overlap, close the circle, rotate so the start
gene (usually dnaA) is base 1 on the forward strand.
"""
import re
from collections import Counter
from .seqs import revcomp


def suffix_prefix_overlap(a, b, k=31, max_ovl=100000, min_ovl=50, max_mismatch=0.02, allow_full=False):
    """Longest overlap where the end of `a` equals the start of `b`
    (substitution-tolerant, no indels).  0 if none."""
    a, b = a.upper(), b.upper()
    if len(a) < k or len(b) < k:
        return 0
    seed = b[:k]
    lo = max(0, len(a) - max_ovl)
    best = 0
    i = a.find(seed, lo)
    while i != -1:
        ovl = len(a) - i
        if ovl >= min_ovl and ovl <= len(b) and (allow_full or ovl < len(a)):
            mm = sum(x != y for x, y in zip(a[i:], b[:ovl]))
            if mm <= max_mismatch * ovl:
                best = max(best, ovl)
        i = a.find(seed, i + 1)
    return best


def end_overlap(seq, **kw):
    """Duplicated suffix == prefix of one contig (circular contig with
    overlapping ends), or 0."""
    kw.setdefault("min_ovl", 100)
    return suffix_prefix_overlap(seq, seq, **kw)


def oriented_end_seqs(seqs, end_a, end_b, w=100000):
    """Sequences laid out so end_a is on the right of the first and end_b on
    the left of the second -- the layout the join will produce."""
    a = seqs[end_a[0]]
    a = a[-w:] if end_a[1] == "R" else revcomp(a[:w])
    b = seqs[end_b[0]]
    b = b[:w] if end_b[1] == "L" else revcomp(b[-w:])
    return a, b


def refine_gap(seqs, end_a, end_b, gap, max_search=5000):
    """Replace a read-based gap estimate near zero (or an overlap) with the
    exact sequence overlap when one exists."""
    if gap is not None and gap > 500:
        return gap
    a, b = oriented_end_seqs(seqs, end_a, end_b)
    ovl = suffix_prefix_overlap(a, b, max_ovl=max(max_search, -(gap or 0) * 2 + 500))
    return -ovl if ovl else gap


def name_hint(name):
    if re.match(r"^(ptg|utg|ctg|htg)\d+c$", name):
        return "hifiasm-style 'c' suffix"
    if re.search(r"circular[=_]?(true|yes|y)\b", name, re.I):
        return "name says circular"
    return None


def find_start(seq, query, k=15, min_votes=8):
    """Locate a start gene (e.g. dnaA from a related strain) on a circle.

    Returns {"pos", "strand", "votes", "identity_hint"} where pos is the
    0-based coordinate of the gene's first base on its strand, or None.
    Tolerates divergence: votes over many k-mers on a shared diagonal."""
    s = seq.upper()
    L = len(s)
    q = query.upper()
    ext = s + s[:len(q)]           # allow a gene that spans the current origin
    index = {}
    for i in range(len(ext) - k + 1):
        index.setdefault(ext[i:i + k], []).append(i)
    best = None
    for strand, qq in (("+", q), ("-", revcomp(q))):
        votes = Counter()
        for i in range(len(qq) - k + 1):
            for j in index.get(qq[i:i + k], ())[:4]:
                votes[(j - i) % L] += 1
        if not votes:
            continue
        diag, n = votes.most_common(1)[0]
        # pool near-diagonals (indels between strains)
        n = sum(votes.get((diag + d) % L, 0) for d in range(-3, 4))
        if n >= min_votes and (best is None or n > best["votes"]):
            pos = diag if strand == "+" else (diag + len(q) - 1) % L
            best = {"pos": pos, "strand": strand, "votes": n,
                    "kmers_in_query": len(qq) - k + 1}
    return best


def circularity_report(graph, candidates, seqs):
    """One row per contig with any circularity signal."""
    rows = []
    links = graph.link_keys()
    by_id = {c["id"]: c for c in candidates}
    for name, seq in seqs.items():
        key_id = f"{name}:L~{name}:R"
        c = by_id.get(key_id)
        self_link = ((name, "L"), (name, "R")) in links
        ovl = end_overlap(seq)
        hint = name_hint(name)
        if not (c or self_link or ovl or hint):
            continue
        rows.append({
            "contig": name, "length": len(seq),
            "graph_self_link": self_link,
            "name_hint": hint,
            "end_overlap": ovl,
            "wrap_longreads": c["longread_reads"] if c else 0,
            "wrap_softclips": c["softclip_reads"] if c else 0,
            "wrap_pairs": c.get("pair_reads", 0) if c else 0,
            "score": c["score"] if c else 0.0,
            "suggested_gap": -ovl if ovl else (c["gap"] if c else 0),
        })
    rows.sort(key=lambda r: -r["score"])
    return rows
