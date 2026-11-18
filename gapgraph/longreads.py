"""Long-read bridges: reads whose alignment hops from one contig end to another.

Consecutive hits along a read (sorted by read coordinate) define a join of
two contig ends plus a gap estimate.  Reads touching >=3 contigs also give
*paths*, which resolve collapsed repeats in the graph (A->R->B vs A->R->C).
"""
from collections import Counter, defaultdict
from .gfa import join_key


def _exit_end(h, slack):
    """End through which the read leaves this hit, plus ref bases left unaligned."""
    if h.strand == "+":
        return (h.ref, "R"), h.rlen - h.re
    return (h.ref, "L"), h.rs


def _entry_end(h, slack):
    if h.strand == "+":
        return (h.ref, "L"), h.rs
    return (h.ref, "R"), h.rlen - h.re


def bridges(hits, min_mapq=5, min_aln=500, end_slack=1000):
    by_read = defaultdict(list)
    for h in hits:
        if h.mapq >= min_mapq and (h.qe - h.qs) >= min(min_aln, 0.8 * h.rlen):
            by_read[h.qname].append(h)
    joins = defaultdict(lambda: {"reads": 0, "gaps": [], "platforms": Counter(), "conflict_reads": 0})
    paths = Counter()
    for q, hs in by_read.items():
        if len(hs) < 2:
            continue
        hs.sort(key=lambda h: h.qs)
        path = []
        for h1, h2 in zip(hs, hs[1:]):
            e1, un1 = _exit_end(h1, end_slack)
            e2, un2 = _entry_end(h2, end_slack)
            if h1.ref == h2.ref and h1.strand == h2.strand and abs(h1.qe - h2.qs) < end_slack:
                wraps = e1 != e2 and un1 <= end_slack and un2 <= end_slack
                if not wraps:
                    continue  # split alignment within one contig (indel/SV), not a bridge
                # else: read runs off one end and continues at the other ->
                # circular molecule wrapping the origin
            key = join_key(e1, e2)
            j = joins[key]
            if un1 > end_slack or un2 > end_slack:
                j["conflict_reads"] += 1  # read leaves mid-contig: possible misassembly
                continue
            j["reads"] += 1
            j["gaps"].append((h2.qs - h1.qe) - un1 - un2)
            j["platforms"][h1.platform] += 1
            if not path:
                path.append((h1.ref, h1.strand))
            path.append((h2.ref, h2.strand))
        if len(path) >= 3:
            # canonicalise direction
            rev = tuple((n, "-" if o == "+" else "+") for n, o in reversed(path))
            paths[min(tuple(path), rev)] += 1
    out = {}
    for key, j in joins.items():
        g = sorted(j["gaps"])
        out[key] = {"reads": j["reads"], "gap": g[len(g) // 2] if g else None,
                    "conflict_reads": j["conflict_reads"], "platforms": dict(j["platforms"]), "ends": key}
    return out, [{"path": [n + o for n, o in p], "reads": n} for p, n in paths.most_common()]
