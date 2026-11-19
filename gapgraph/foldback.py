"""Foldbacks: reads and contig ends that run into their own reverse complement.

Foldback READS (mostly nanopore: the complementary strand of the same
molecule follows the template through the pore without a separating read
break) carry a sequence followed by its reverse complement, usually with a
short hairpin/adapter between the arms.  Aligners report them as
  * a primary + supplementary pair on opposite strands of the same place, or
  * one alignment whose soft clip is the reverse complement of the
    reference just before the clip point.
Left in, the folded arm double-counts depth and masquerades as soft-clip /
long-read evidence for inversions, mis-joins and false joins.

Foldback CONTIG ENDS (assembler hairpins: a contig whose last n bases are the
reverse complement of the n before them) are artefacts to trim -- unless
reads genuinely run straight through the fold, as on linear replicons with
hairpin telomeres (e.g. Borrelia linear plasmids).  So the read evidence is
reported with every call, and removal is a recommendation, not automatic,
when reads span the fold.

Modes for foldback reads:
    keep   leave alignments as they are (for inspection)
    trim   keep the longer arm, hard-clip the folded arm (default)
    drop   remove every alignment of the read
"""
from collections import Counter, defaultdict
from dataclasses import replace

from .alignments import AlignmentStore
from .realign import prefix_align
from .seqs import revcomp


def _key(a):
    return (a.qname, a.ref, a.pos, a.reverse, a.supplementary)


def _aligned(a):
    return sum(n for n, op in a.cigar if op in "M=X")


def _kmer_containment(query, target, k=15):
    q = query.upper()
    t = target.upper()
    if len(q) < k or len(t) < k:
        return 0.0
    ts = {t[i:i + k] for i in range(len(t) - k + 1)}
    n = len(q) - k + 1
    return sum(1 for i in range(n) if q[i:i + k] in ts) / n


# ------------------------------------------------------------------ reads
def detect_read_foldbacks(store, seqs, max_ref_gap=500, max_query_gap=1000, min_overlap_frac=0.3,
                          min_clip=200, min_contain=0.15, exclude=None, pile=3, pile_window=50, stats=None):
    """Return a list of foldback events, one per read.

    Reads that only *look* folded because the reference is (inside the arm of
    a contig foldback, `exclude` = {contig: [(lo, hi)]}), and piles of >= `pile`
    folded reads at one position (a palindrome in the reference: real read
    foldbacks fall at random positions) are not counted; their numbers go
    into `stats` if a dict is given."""
    by_q = defaultdict(list)
    for a in store.all():
        by_q[a.qname].append(a)
    events = []
    for q, alns in by_q.items():
        ev = None
        if len(alns) >= 2:
            # opposite-strand arms of one read on the same stretch of one contig
            for i in range(len(alns)):
                for j in range(i + 1, len(alns)):
                    a, b = alns[i], alns[j]
                    if a.ref != b.ref or a.reverse == b.reverse or (a.paired and a.flag & 192 != b.flag & 192):
                        continue
                    ov = min(a.end, b.end) - max(a.pos, b.pos)
                    short = min(a.end - a.pos, b.end - b.pos)
                    if not (ov >= min_overlap_frac * short or -ov <= max_ref_gap):
                        continue
                    (qa0, qa1), (qb0, qb1) = a.query_interval(), b.query_interval()
                    first, second = ((qa0, qa1), (qb0, qb1)) if qa0 <= qb0 else ((qb0, qb1), (qa0, qa1))
                    qgap = second[0] - first[1]
                    if not (-50 <= qgap <= max_query_gap):
                        continue
                    keep, arm = (a, b) if _aligned(a) >= _aligned(b) else (b, a)
                    k0, k1 = keep.query_interval()
                    a0, _ = arm.query_interval()
                    three_prime = a0 >= k1 - 50          # folded arm follows the kept arm in the read
                    # read 3' end is on the right for a forward alignment
                    side = ("R" if not keep.reverse else "L") if three_prime else ("L" if not keep.reverse else "R")
                    ev = {"qname": q, "kind": "split", "ref": keep.ref, "side": side,
                          "pos": keep.end if side == "R" else keep.pos,
                          "arm": arm.end - arm.pos, "loop": qgap, "keep": _key(keep), "drop": [_key(arm)],
                          "platform": keep.platform}
                    break
                if ev:
                    break
        if ev is None:
            # a single alignment whose clip is the reverse complement of what precedes it
            for a in alns:
                if not a.seq or a.ref not in seqs:
                    continue
                ref = seqs[a.ref]
                if a.rclip_soft >= min_clip:
                    clip = a.seq[len(a.seq) - a.rclip_soft:]
                    region = ref[max(0, a.end - len(clip) - 200):a.end]
                    if _kmer_containment(clip, revcomp(region)) >= min_contain:
                        ev = {"qname": q, "kind": "clip", "ref": a.ref, "side": "R", "pos": a.end,
                              "arm": a.rclip_soft, "keep": _key(a), "drop": [], "platform": a.platform}
                        break
                if a.lclip_soft >= min_clip:
                    clip = a.seq[:a.lclip_soft]
                    region = ref[a.pos:a.pos + len(clip) + 200]
                    if _kmer_containment(clip, revcomp(region)) >= min_contain:
                        ev = {"qname": q, "kind": "clip", "ref": a.ref, "side": "L", "pos": a.pos,
                              "arm": a.lclip_soft, "keep": _key(a), "drop": [], "platform": a.platform}
                        break
        if ev:
            events.append(ev)
    exclude = exclude or {}
    kept, at_contig = [], 0
    for e in events:
        if any(lo <= e["pos"] <= hi for lo, hi in exclude.get(e["ref"], ())):
            at_contig += 1
        else:
            kept.append(e)
    piles = Counter((e["ref"], e["pos"] // pile_window) for e in kept)
    near = lambda e: sum(piles.get((e["ref"], e["pos"] // pile_window + d), 0) for d in (-1, 0, 1))
    final = [e for e in kept if near(e) < pile]
    if stats is not None:
        stats.update(at_contig_foldback=at_contig, reference_palindrome_piles=len(kept) - len(final))
    return final


def hard_clip(a, side):
    """Turn the soft clip on one side into a hard clip (bases leave the read)."""
    cig = list(a.cigar)
    if side == "R" and a.rclip_soft:
        n = a.rclip_soft
        h = cig.pop()[0] if cig[-1][1] == "H" else 0
        cig.pop()                                  # the S op
        cig.append((n + h, "H"))
        return replace(a, seq=a.seq[:-n], qual=a.qual[:-n] if a.qual else "", cigar=cig, rclip_soft=0)
    if side == "L" and a.lclip_soft:
        n = a.lclip_soft
        h = cig.pop(0)[0] if cig[0][1] == "H" else 0
        cig.pop(0)
        cig.insert(0, (n + h, "H"))
        return replace(a, seq=a.seq[n:], qual=a.qual[n:] if a.qual else "", cigar=cig, lclip_soft=0)
    return a


def apply_mode(store, events, mode="trim"):
    """A new AlignmentStore with foldback reads handled per `mode`."""
    if mode == "keep" or not events:
        return store
    by_q = {e["qname"]: e for e in events}
    drop = {k for e in events for k in e["drop"]}
    out = []
    for a in store.all():
        e = by_q.get(a.qname)
        if e is None:
            out.append(a)
            continue
        if mode == "drop":
            continue
        k = _key(a)
        if k in drop:
            continue
        if k == e["keep"]:
            a = hard_clip(a, e["side"])
            if a.supplementary:     # the kept arm becomes the primary
                a = replace(a, supplementary=False, flag=a.flag & ~2048)
        out.append(a)
    return AlignmentStore().add(out)


def summarise(events, store, stats=None):
    n_reads = len({a.qname for a in store.all()})
    by = Counter(e["kind"] for e in events)
    return {**(stats or {}), "reads": n_reads, "foldbacks": len(events), "fraction": round(len(events) / max(1, n_reads), 4),
            "split": by.get("split", 0), "clip": by.get("clip", 0),
            "median_arm": sorted(e["arm"] for e in events)[len(events) // 2] if events else 0}


# ------------------------------------------------------------ contig ends
def end_foldback(seq, side, window=60000, k=21, min_arm=300, min_kmers=15, min_ident=0.85):
    """Hairpin at one contig end: the last n bases are the reverse complement
    of the n before them (possibly with a short loop between).

    Works on the end as the right-hand side of `s`.  A fold at f in window
    S makes S[f + t] == comp(S[f - 1 - t]), i.e. S[i] == R[i + d] with
    R = revcomp(S) and d = W - 2f: one diagonal of S-vs-R k-mer matches that
    runs to the end of the window."""
    s = seq if side == "R" else revcomp(seq)
    L = len(s)
    W = min(L, window)
    S = s[L - W:].upper()
    R = revcomp(S)
    idx = defaultdict(list)
    for j in range(W - k + 1):
        km = R[j:j + k]
        if len(idx[km]) < 4:
            idx[km].append(j)
    diag = Counter()
    last = {}
    for i in range(W - k + 1):
        for j in idx.get(S[i:i + k], ()):
            d = j - i
            diag[d] += 1
            if i > last.get(d, -1):
                last[d] = i
    best = None
    for d, n in diag.most_common(20):
        if n < min_kmers or last[d] < W - k - 300:
            continue
        f = (W - d) // 2
        arm = W - f
        if arm < min_arm or f <= 0:
            continue
        al = prefix_align(S[f:], revcomp(S[:f]))
        if al["identity"] >= min_ident and al["q"] >= 0.8 * arm:
            cand = {"fold_s": f, "arm": arm, "identity": round(al["identity"], 4), "kmers": n}
            if best is None or cand["arm"] > best["arm"]:
                best = cand
    if best is None:
        return None
    p = L - W + best.pop("fold_s")           # fold point in s coordinates
    fold = p if side == "R" else L - p       # keep [0, fold) for R; keep [fold, L) for L
    return {"side": side, "fold": fold, **best}


def contig_foldbacks(seqs, stores, skip=(), events=(), hits=None):
    """Scan every contig end; attach read evidence at each fold point.
    `hits` (precomputed end_foldback results) skips the sequence scan."""
    fb_reads = {e["qname"] for e in events}
    out = []
    for name, seq in seqs.items():
        if name in skip:
            continue
        for hit in (hits if hits is not None else (end_foldback(seq, s) for s in "RL")):
            if not hit:
                continue
            x = hit["fold"]
            span = depth = 0
            for st in stores.values():
                for a in st.fetch(name, max(0, x - 1), x + 1):
                    if a.supplementary or a.qname in fb_reads:
                        continue
                    depth += 1
                    if a.pos < x - 50 and a.end > x + 50:
                        span += 1
            hit.update(contig=name, length=len(seq), spanning_reads=span, depth=depth)
            # an assembler hairpin has no genuine reads through the fold; a real
            # hairpin telomere / terminal inverted repeat does
            hit["recommend"] = "remove" if span <= max(2, 0.1 * depth) else "check"
            out.append(hit)
    return out
