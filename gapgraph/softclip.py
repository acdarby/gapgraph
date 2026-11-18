"""Soft-clip evidence: what Gap5 users did by eye, done systematically.

1. At contig *ends*, reads that run off the end carry soft-clipped sequence
   that belongs to whatever comes next.  We orient every clip "outward" and
   k-mer match it against the *inward* sequence of every other contig end
   (both orientations).  A consistent k-mer diagonal gives partner + gap.
2. *Inside* contigs, piles of reads clipping at the same column are the
   classic signature of a mis-join / collapsed repeat boundary -> break
   candidates.
"""
from collections import Counter, defaultdict
from .seqs import revcomp, kmer_index


def inward_seq(seq, side, w):
    """Sequence read from the end *into* the contig."""
    return seq[:w] if side == "L" else revcomp(seq[-w:])


def outward_clips(store, contig_len, min_clip, end_window, min_mapq):
    """{(contig, side): [(qname, outward_clip_seq, platform)]}"""
    out = defaultdict(list)
    for a in store.all():
        if a.mapq < min_mapq:
            continue
        L = contig_len.get(a.ref)
        if L is None:
            continue
        if a.lclip_soft >= min_clip and a.pos <= end_window:
            # clip precedes the contig's L end; outward = continue on - strand
            out[(a.ref, "L")].append((a.qname, revcomp(a.lclip_seq[:a.lclip_soft]), a.platform))
        if a.rclip_soft >= min_clip and a.end >= L - end_window:
            out[(a.ref, "R")].append((a.qname, a.rclip_seq, a.platform))
    return out


def internal_hotspots(store, contig_len, min_clip, end_window, min_mapq, min_reads=3, bin_size=20):
    hits = defaultdict(Counter)
    exact = defaultdict(list)
    depth = defaultdict(int)
    for a in store.all():
        if a.mapq < min_mapq:
            continue
        L = contig_len.get(a.ref, 0)
        depth[a.ref] += a.end - a.pos
        if a.lclip >= min_clip and a.pos > end_window:
            hits[a.ref][(a.pos // bin_size, "left")] += 1
            exact[(a.ref, a.pos // bin_size, "left")].append(a.pos)
        if a.rclip >= min_clip and a.end < L - end_window:
            hits[a.ref][(a.end // bin_size, "right")] += 1
            exact[(a.ref, a.end // bin_size, "right")].append(a.end)
    res = []
    for ref, c in hits.items():
        mean_cov = depth[ref] / max(1, contig_len.get(ref, 1))
        thr = max(min_reads, 0.25 * mean_cov)
        for (b, side), n in c.items():
            if n >= thr:
                xs = sorted(exact[(ref, b, side)])
                res.append({"contig": ref, "pos": xs[len(xs) // 2], "clip_side": side,
                            "reads": n, "mean_cov": round(mean_cov, 1)})
    res.sort(key=lambda h: -h["reads"])
    return res


def match_clips(clips, seqs, k, max_inward=50000, min_shared=4, self_ok=False):
    """Return {join_key: {"reads": n, "gap": median_gap, "ends": (a, b)}}."""
    ends = [(n, s) for n in seqs for s in "LR"]
    index = defaultdict(list)  # kmer -> [(end, pos)]
    for e in ends:
        inw = inward_seq(seqs[e[0]], e[1], max_inward).upper()
        for km, ps in kmer_index(inw, k).items():
            if len(ps) <= 2:
                for p in ps:
                    index[km].append((e, p))
    support = defaultdict(list)
    for end_a, lst in clips.items():
        for qname, clip, plat in lst:
            c = clip.upper()
            diag = Counter()
            for i in range(0, len(c) - k + 1):
                for e, p in index.get(c[i:i + k], ()):
                    if e == end_a or (not self_ok and e[0] == end_a[0] and e[1] == end_a[1]):
                        continue
                    diag[(e, (i - p) // 25)] += 1   # gap bucket of 25bp tolerates indels
            if not diag:
                continue
            (e, gb), n = diag.most_common(1)[0]
            if n >= min_shared:
                support[(end_a, e)].append((gb * 25, qname, plat))
    from .gfa import join_key
    merged = {}
    for (a, b), lst in support.items():
        key = join_key(a, b)
        m = merged.setdefault(key, {"reads": 0, "gaps": [], "platforms": Counter(), "ends": key})
        m["reads"] += len(lst)
        m["gaps"].extend(g for g, _, _ in lst)
        m["platforms"].update(p for _, _, p in lst)
    for m in merged.values():
        g = sorted(m.pop("gaps"))
        m["gap"] = g[len(g) // 2]
        m["platforms"] = dict(m["platforms"])
    return merged


def naive_consensus(clip_seqs, max_len=2000, min_depth=2):
    """Column-majority consensus of outward clips (good for HiFi/Illumina,
    rough for ONT/CLR where indels misalign columns)."""
    cons = []
    for i in range(max_len):
        col = Counter(s[i].upper() for s in clip_seqs if len(s) > i)
        if sum(col.values()) < min_depth:
            break
        cons.append(col.most_common(1)[0][0])
    return "".join(cons)
