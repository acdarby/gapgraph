"""Read pairs and templates (Gap5 "Find Read Pairs" + "Template Display").

Template status, using Gap5's colour scheme:
    consistent    FR pair, insert within the library's normal range (grey, by MAPQ)
    size          FR pair, insert too short or too long                (grey-red)
    orientation   both reads in the same orientation, or pointing away (red)
    spanning      mate (or, for long reads, another alignment) in another contig (orange)
    wrap          outward-facing pair across the two ends of one contig: circular (green)
    single        no mate in the data (blue)

Pairs whose mates sit in different contigs near the contig ends are
scaffolding evidence: they give the two joined ends and a gap estimate.
"""
from collections import Counter, defaultdict
import random

from .gfa import join_key

MAD_K = 1.4826


def insert_stats(store, min_mapq=20, max_n=200000):
    """Median and robust SD of insert size for FR pairs, per platform."""
    sizes = defaultdict(list)
    for a in store.all():
        if (a.paired and not a.supplementary and a.mapq >= min_mapq and a.mate_ref == a.ref
                and a.tlen > 0 and not a.reverse and a.mate_reverse and a.mate_pos >= a.pos):
            if len(sizes[a.platform]) < max_n:
                sizes[a.platform].append(a.tlen)
    out = {}
    for plat, v in sizes.items():
        if len(v) < 20:
            continue
        v.sort()
        med = v[len(v) // 2]
        mad = sorted(abs(x - med) for x in v)[len(v) // 2] * MAD_K
        out[plat] = {"median": med, "sd": max(mad, 0.05 * med), "n": len(v),
                     "lo": max(0, med - 4 * max(mad, 0.05 * med)), "hi": med + 4 * max(mad, 0.05 * med)}
    return out


def classify(a, stats, other_contigs=None):
    """Status of the template this alignment belongs to."""
    if a.paired and a.mate_ref:
        if a.mate_ref != a.ref:
            return "spanning"
        st0 = stats.get(a.platform)
        if st0 and abs(a.tlen) > st0["hi"]:
            outward = (not a.reverse and a.mate_reverse and a.mate_pos < a.pos) or \
                      (a.reverse and not a.mate_reverse and a.mate_pos > a.pos)
            if outward:
                return "wrap"   # spans the end of a circular contig
        st = stats.get(a.platform)
        left_fwd = (not a.reverse) if a.pos <= a.mate_pos else (not a.mate_reverse)
        right_rev = a.mate_reverse if a.pos <= a.mate_pos else a.reverse
        if a.reverse == a.mate_reverse or not (left_fwd and right_rev):
            return "orientation"
        if st and not (st["lo"] <= abs(a.tlen) <= st["hi"]):
            return "size"
        return "consistent"
    if other_contigs and other_contigs.get(a.qname, set()) - {a.ref}:
        return "spanning"
    return "single"


def spanning_joins(store, contig_len, stats, min_mapq=20):
    """Pairs with mates on different contigs -> end joins with gap estimates."""
    acc = defaultdict(lambda: {"pairs": 0, "gaps": []})
    for a in store.all():
        if not (a.paired and a.mate_ref and a.first_in_pair and not a.supplementary and a.mapq >= min_mapq):
            continue
        if a.mate_ref == a.ref:
            # same contig: only interesting if the pair points *outward* across
            # the two ends -- the signature of a circular molecule
            wrap = (not a.reverse and a.mate_reverse and a.mate_pos < a.pos) or \
                   (a.reverse and not a.mate_reverse and a.mate_pos > a.pos)
            if not wrap:
                continue
        st = stats.get(a.platform)
        if not st:
            continue
        lx, ly = contig_len.get(a.ref), contig_len.get(a.mate_ref)
        if lx is None or ly is None:
            continue
        rl = a.end - a.pos
        # the read points out of its contig through this end:
        if not a.reverse:
            ex, dx = (a.ref, "R"), lx - a.pos
        else:
            ex, dx = (a.ref, "L"), a.end
        # mate pointing back towards the read enters through:
        if a.mate_reverse:
            ey, dy = (a.mate_ref, "L"), a.mate_pos + rl
        else:
            ey, dy = (a.mate_ref, "R"), ly - a.mate_pos
        if dx > st["hi"] or dy > st["hi"]:
            continue  # too far from the ends to bridge
        j = acc[join_key(ex, ey)]
        j["pairs"] += 1
        j["gaps"].append(int(st["median"] - dx - dy))
    out = {}
    for k, j in acc.items():
        g = sorted(j["gaps"])
        out[k] = {"pairs": j["pairs"], "gap": g[len(g) // 2]}
    return out


def other_contig_index(stores):
    """qname -> set of contigs it aligns to (for long-read 'spanning')."""
    idx = defaultdict(set)
    for st in stores.values():
        for a in st.all():
            if not a.paired:
                idx[a.qname].add(a.ref)
    return idx


def templates(stores, contig, start, end, stats, other, max_items=6000, seed=1):
    """Template-display items in [start, end)."""
    items = []
    for plat, st in stores.items():
        for a in st.fetch(contig, start, end):
            if a.supplementary and a.paired:
                continue
            status = classify(a, stats, other)
            if a.paired and a.mate_ref == a.ref and status != "orientation":
                if a.pos > a.mate_pos or (a.pos == a.mate_pos and not a.first_in_pair):
                    continue  # draw each pair once, from its leftmost read
                x0, x1 = a.pos, max(a.end, a.mate_pos + (a.end - a.pos))
                items.append({"x0": x0, "x1": x1, "y": x1 - x0, "s": status, "q": a.mapq, "p": plat,
                              "r": [[a.pos, a.end, 0], [a.mate_pos, a.mate_pos + a.end - a.pos, 1]], "n": a.qname})
            else:
                items.append({"x0": a.pos, "x1": a.end, "y": a.end - a.pos, "s": status, "q": a.mapq,
                              "p": plat, "r": [[a.pos, a.end, 1 if a.reverse else 0]], "n": a.qname})
    counts = Counter(i["s"] for i in items)
    if len(items) > max_items:
        rng = random.Random(seed)
        keep = [i for i in items if i["s"] != "consistent" and i["s"] != "single"]
        rest = [i for i in items if i["s"] in ("consistent", "single")]
        k = max(0, max_items - len(keep))
        items = keep[:max_items] + rng.sample(rest, min(k, len(rest)))
    return {"items": items, "counts": dict(counts), "sampled": sum(counts.values()) > len(items)}


def depth_track(stores, contig, length, start, end, bins=600):
    """Mean read depth and template (pair-span) depth per bin."""
    import numpy as np
    start, end = max(0, start), min(length, end)
    n = max(1, end - start)
    rd = np.zeros(n + 1)
    td = np.zeros(n + 1)
    for st in stores.values():
        for a in st.fetch(contig, start, end):
            if a.supplementary:
                continue
            rd[max(0, a.pos - start)] += 1
            rd[min(n, a.end - start)] -= 1
            if a.paired and a.mate_ref == a.ref and a.pos <= a.mate_pos and a.tlen > 0:
                td[max(0, a.pos - start)] += 1
                td[min(n, a.mate_pos + (a.end - a.pos) - start)] -= 1
    rd, td = np.cumsum(rd)[:n], np.cumsum(td)[:n]
    step = max(1, n // bins)
    k = n // step
    r = rd[:k * step].reshape(k, step).mean(axis=1)
    t = td[:k * step].reshape(k, step).mean(axis=1)
    return {"start": start, "step": step, "reads": r.round(1).tolist(), "templates": t.round(1).tolist()}
