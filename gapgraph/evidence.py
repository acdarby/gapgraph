"""Merge all evidence onto join candidates (unordered pairs of contig ends)."""
import math
from collections import defaultdict


def _via_from_paths(graph, lr_paths):
    """Long-read paths X -> R -> Y through a tangle segment R give an X~Y join
    that should carry a copy of R (repeat resolution)."""
    from .gfa import join_key, end_from_orient
    tangles = set(graph.tangles())
    via = {}
    for p in lr_paths:
        steps = [(s[:-1], s[-1]) for s in p["path"]]
        for x, r, y in zip(steps, steps[1:], steps[2:]):
            if r[0] not in tangles or x[0] in tangles or y[0] in tangles:
                continue
            ea, eb = end_from_orient(*x, True), end_from_orient(*y, False)
            key = join_key(ea, eb)
            mid = [r]
            if key[0] != ea:  # key stored b->a: walk the path backwards
                mid = [(r[0], "-" if r[1] == "+" else "+")]
            v = via.setdefault(key, {"via": [n + o for n, o in mid], "reads": 0})
            v["reads"] += p["reads"]
    return via


def build_candidates(graph, clip_joins, lr_joins, hic=None, lr_paths=()):
    links = graph.link_keys()
    via = _via_from_paths(graph, lr_paths)
    keys = set(links) | set(clip_joins) | set(lr_joins) | set(via)
    if hic is not None:
        keys |= hic.candidate_keys()
    cands = []
    for key in keys:
        (a, sa), (b, sb) = key
        if a == b and sa == sb:
            continue
        link = links.get(key)
        sc = clip_joins.get(key, {})
        lr = lr_joins.get(key, {})
        h = hic.score(key) if hic is not None else None
        c = {
            "id": f"{a}:{sa}~{b}:{sb}",
            "a": [a, sa], "b": [b, sb],
            "graph_link": link is not None,
            "graph_overlap": link.overlap if link else None,
            "softclip_reads": sc.get("reads", 0),
            "softclip_gap": sc.get("gap"),
            "softclip_platforms": sc.get("platforms", {}),
            "longread_reads": lr.get("reads", 0),
            "longread_gap": lr.get("gap"),
            "longread_platforms": lr.get("platforms", {}),
            "hic_contacts": h["contacts"] if h else 0,
            "hic_norm": h["norm"] if h else 0.0,
            "hic_orient": h["orient_frac"] if h else 0.0,
            "flags": [],
            "via": via.get(key, {}).get("via", []),
            "via_reads": via.get(key, {}).get("reads", 0),
        }
        s = 0.0
        s += 3.0 if c["graph_link"] else 0.0
        s += 2.0 * math.log2(1 + max(c["longread_reads"], c["via_reads"]))
        s += 1.0 * math.log2(1 + c["softclip_reads"])
        if h and h["contacts"]:
            s += 6.0 * c["hic_norm"] * (0.5 + c["hic_orient"])
        c["score"] = round(s, 2)
        gaps = [g for g in (c["longread_gap"], c["softclip_gap"]) if g is not None]
        c["gap"] = gaps[0] if gaps else (-c["graph_overlap"] if c["graph_overlap"] else 100)
        if c["via"]:
            c["gap"] = 0
        c["evidence"] = [n for n, v in (("graph", c["graph_link"]), ("longread", c["longread_reads"]),
                                        ("softclip", c["softclip_reads"]), ("hic", c["hic_contacts"] >= 5),
                                        ("path", c["via_reads"])) if v]
        cands.append(c)
    _flag(cands)
    cands.sort(key=lambda c: -c["score"])
    return cands


def _partner(c, end):
    return tuple(c["b"]) if tuple(c["a"]) == end else tuple(c["a"])


def _consistent(c1, c2, end):
    """X~R and X~Y(via R) agree: the second just continues through the repeat."""
    for x, y in ((c1, c2), (c2, c1)):
        if x["via"] and _partner(y, end)[0] in [v[:-1] for v in x["via"]]:
            return True
    return False


def _flag(cands):
    per_end = defaultdict(list)
    for c in cands:
        per_end[tuple(c["a"])].append(c)
        per_end[tuple(c["b"])].append(c)
    for end, lst in per_end.items():
        lst.sort(key=lambda c: -c["score"])
        best = lst[0]
        for c in lst[1:]:
            if _consistent(best, c, end):
                continue
            if c["score"] >= 0.6 * best["score"]:
                for x in (best, c):
                    if "ambiguous" not in x["flags"]:
                        x["flags"].append("ambiguous")
            elif best is not c and c["graph_link"] and "outcompeted" not in c["flags"]:
                c["flags"].append("outcompeted")  # graph says yes, reads/Hi-C prefer another partner
    for c in cands:
        if c["via"]:
            c["flags"].append("via:" + ",".join(c["via"]))
        if c["evidence"] == ["hic"]:
            c["flags"].append("hic-only")
        if c["longread_gap"] is not None and c["longread_gap"] > 0 and not c["graph_link"]:
            c["flags"].append("gap-to-fill")
        if c["graph_link"] and not (c["longread_reads"] or c["softclip_reads"]) and c["hic_contacts"] < 5:
            c["flags"].append("unsupported-link")
