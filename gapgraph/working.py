"""The working assembly: what the editor shows after every edit.

Gap5 edits contigs in place: Break Contig leaves two contigs with the reads
divided between them, and a join leaves one contig with all the reads.
gapgraph keeps the input alignments read-only, so after every operation it
*derives* the working assembly from the project layout:

  * each scaffold becomes one working contig (sequence = its components,
    overlaps trimmed, circles rotated)
  * every input alignment is cut to the component(s) it falls in (bases
    outside become soft clips, exactly what an aligner would report against
    the new contig) and lifted into working coordinates, reverse-complemented
    for components used on the minus strand; mates are lifted too
  * at internal junctions created by joins, the soft clips that ran off the
    old contig end are re-aligned into the sequence on the other side, so
    reads now span the join, as they would in a joined Gap5 contig
  * graph links, PAF hits and Hi-C contacts are mapped to working contigs

Scaffolds whose single part is an unflipped whole input contig reuse the
input alignments untouched, so an unedited project costs nothing.
"""
from bisect import bisect_right
from dataclasses import replace

from .alignments import Alignment, Hit
from .gfa import AssemblyGraph, Link, Segment
from .realign import prefix_align, cols_to_cigar
from .seqs import revcomp


def working_name(s):
    parts = s["parts"]
    if len(parts) == 1 and not s.get("circular"):
        return parts[0]["piece"]
    return s["name"]


# --------------------------------------------------------------- alignment ops
def _expand(a):
    """Per-column events: (op, ref_pos) with op in M, I, D."""
    out = []
    r = a.pos
    for n, op in a.cigar:
        if op in "M=X":
            out.extend(("M", r + i) for i in range(n))
            r += n
        elif op == "I":
            out.extend(("I", None) for _ in range(n))
        elif op in "DN":
            out.extend(("D", r + i) for i in range(n))
            r += n
    return out


def _compress(events):
    cig = []
    for op, _ in events:
        if cig and cig[-1][1] == op:
            cig[-1][0] += 1
        else:
            cig.append([1, op])
    return [(n, op) for n, op in cig]


def _lead_hard(a):
    return a.cigar[0][0] if a.cigar and a.cigar[0][1] == "H" else 0


def _trail_hard(a):
    return a.cigar[-1][0] if a.cigar and a.cigar[-1][1] == "H" else 0


def _rebuild(a, pos, core, lsoft, rsoft):
    lh, th = _lead_hard(a), _trail_hard(a)
    cig = ([(lh, "H")] if lh else []) + ([(lsoft, "S")] if lsoft else []) + core + \
          ([(rsoft, "S")] if rsoft else []) + ([(th, "H")] if th else [])
    rlen = sum(n for n, op in core if op in "MD=XN")
    return replace(a, pos=pos, end=pos + rlen, cigar=cig, lclip=lh + lsoft, rclip=th + rsoft,
                   lclip_soft=lsoft, rclip_soft=rsoft)


def cut_alignment(a, lo, hi):
    """Restrict an alignment to reference [lo, hi); bases outside become soft clips."""
    if a.pos >= lo and a.end <= hi:
        return a
    if a.end <= lo or a.pos >= hi or not a.seq:
        return None
    ev = _expand(a)
    keep = [i for i, (op, r) in enumerate(ev) if op == "M" and lo <= r < hi]
    if not keep:
        return None
    f, l = keep[0], keep[-1]
    q_before = sum(1 for op, _ in ev[:f] if op != "D")
    q_after = sum(1 for op, _ in ev[l + 1:] if op != "D")
    core = _compress(ev[f:l + 1])
    return _rebuild(a, ev[f][1], core, a.lclip_soft + q_before, a.rclip_soft + q_after)


def reverse_alignment(a, L):
    """The same alignment against the reverse complement of a length-L reference."""
    core = [(n, op) for n, op in a.cigar if op not in "SH"][::-1]
    lh, th = _lead_hard(a), _trail_hard(a)
    b = replace(a, seq=revcomp(a.seq) if a.seq else "", qual=a.qual[::-1] if a.qual else "",
                reverse=not a.reverse, flag=a.flag ^ 16)
    b = _rebuild(replace(b, cigar=([(th, "H")] if th else []) + [(1, "M")] + ([(lh, "H")] if lh else [])),
                 L - a.end, core, a.rclip_soft, a.lclip_soft)
    return b


# ------------------------------------------------------------------ the build
class Mapper:
    """input contig position -> working (contig, position, strand)."""

    def __init__(self):
        self.by_src = {}

    def add(self, src, a, b, orient, off, wname):
        self.by_src.setdefault(src, []).append((a, b, orient, off, wname))

    def finish(self):
        for v in self.by_src.values():
            v.sort()
        self._starts = {k: [x[0] for x in v] for k, v in self.by_src.items()}

    def map(self, src, p):
        v = self.by_src.get(src)
        if not v:
            return None
        i = bisect_right(self._starts[src], p) - 1
        # several components can share an input region (repeat copies);
        # scan back for the first one that contains p
        while i >= 0:
            a, b, o, off, w = v[i]
            if a <= p < b:
                return (w, off + (p - a), "+") if o == "+" else (w, off + (b - 1 - p), "-")
            i -= 1
        return None


class Working:
    def __init__(self, project, orig_graph, orig_seqs, orig_stores, orig_hits=(), hic_pairs=(),
                 cache=None, realign_junctions=True):
        self.project = project
        self.cache = cache if cache is not None else {}
        self.seqs, self.comps, self.circular, self.scaffold_of = {}, {}, set(), {}
        self.mapper = Mapper()
        for s in project.scaffolds:
            w = working_name(s)
            comps = project.components(s, rotate=True)
            off, placed = 0, []
            for c in comps:
                n = c[1] if c[0] == "N" else c[3] - c[2]
                if c[0] == "W":
                    self.mapper.add(c[1], c[2], c[3], c[4], off, w)
                    placed.append((off, c))
                off += n
            self.seqs[w] = project.scaffold_seq(s)
            self.comps[w] = placed
            self.scaffold_of[w] = s["name"]
            if s.get("circular"):
                self.circular.add(w)
        self.mapper.finish()
        self.lengths = {w: len(s) for w, s in self.seqs.items()}
        self.graph = self._graph(orig_graph)
        self.stores = self._stores(orig_stores, realign_junctions)
        self.hits = [h for h in (self._map_hit(h) for h in orig_hits) if h]
        self.hic_pairs = [p for p in (self._map_pair(x) for x in hic_pairs) if p]

    # ............................................................. identity
    def _identity(self, w, placed, orig_seqs_len):
        if len(placed) != 1:
            return False
        off, c = placed[0]
        return off == 0 and c[4] == "+" and c[2] == 0 and c[3] == orig_seqs_len.get(c[1], -1) and c[1] == w

    # ................................................................ reads
    def _stores(self, orig_stores, realign_junctions):
        from .alignments import AlignmentStore
        src_len = {n: len(s) for n, s in self.project.seqs.items()}
        out = {}
        for plat, st in orig_stores.items():
            new = AlignmentStore()
            alns = []
            for w, placed in self.comps.items():
                key = (plat, w, tuple((off, tuple(c)) for off, c in placed), len(self.seqs[w]))
                if key in self.cache:
                    alns.extend(self.cache[key])
                    continue
                if self._identity(w, placed, src_len):
                    lifted = list(st.by_ref.get(w, []))
                else:
                    lifted = []
                    for off, c in placed:
                        src, a0, b0, o = c[1], c[2], c[3], c[4]
                        for a in st.fetch(src, a0, b0):
                            x = cut_alignment(a, a0, b0)
                            if x is None:
                                continue
                            if o == "+":
                                x = replace(x, pos=x.pos - a0 + off, end=x.end - a0 + off, ref=w)
                            else:
                                x = reverse_alignment(replace(x, pos=x.pos - a0, end=x.end - a0), b0 - a0)
                                x = replace(x, pos=x.pos + off, end=x.end + off, ref=w)
                            lifted.append(x)
                    if realign_junctions and len(placed) > 1:
                        lifted = self._realign_junctions(w, placed, lifted)
                    lifted = _dedupe(lifted)
                self.cache[key] = lifted
                alns.extend(lifted)
            alns = [self._lift_mate(a) for a in alns]
            out[plat] = new.add(alns)
        return out

    def _lift_mate(self, a):
        if not (a.paired and a.mate_ref):
            return a
        rl = a.end - a.pos
        m = self.mapper.map(a.mate_ref, a.mate_pos)
        if m is None:
            return replace(a, mate_ref="", flag=a.flag | 8)
        w, mp, o = m
        flag = a.flag
        if o == "-":
            mp = mp - rl + 1      # leftmost of the mate on the reversed component
            flag ^= 32
        tlen = 0
        if w == a.ref:
            lo, hi = min(a.pos, mp), max(a.end, mp + rl)
            tlen = (hi - lo) if a.pos <= mp else -(hi - lo)
        return replace(a, mate_ref=w, mate_pos=max(0, mp), tlen=tlen, flag=flag)

    def _realign_junctions(self, w, placed, alns, slack=25, min_ident=0.75):
        """Extend reads clipped at an internal junction into the other side."""
        seq = self.seqs[w]
        junctions = set()
        for (o1, c1), (o2, c2) in zip(placed, placed[1:]):
            e1 = o1 + c1[3] - c1[2]
            if e1 == o2:                      # adjacent, no N gap between
                junctions.add(e1)
        if not junctions:
            return alns
        js = sorted(junctions)
        out = []
        for a in alns:
            if a.seq and a.rclip_soft >= 10:
                k = bisect_right(js, a.end + slack) - 1
                if k >= 0 and js[k] >= a.end - slack:
                    clip = a.seq[len(a.seq) - a.rclip_soft:]
                    tgt = seq[a.end:a.end + len(clip) + 100]
                    r = prefix_align(clip, tgt)
                    if r["t"] >= 10 and r["identity"] >= min_ident:
                        core = [(n, op) for n, op in a.cigar if op not in "SH"] + cols_to_cigar(r["cols"])
                        a = _rebuild(a, a.pos, _merge(core), a.lclip_soft, a.rclip_soft - r["q"])
            if a.seq and a.lclip_soft >= 10:
                k = bisect_right(js, a.pos + slack) - 1
                if k >= 0 and a.pos - slack <= js[k] <= a.pos + slack:
                    clip = revcomp(a.seq[:a.lclip_soft])
                    tgt = revcomp(seq[max(0, a.pos - a.lclip_soft - 100):a.pos])
                    r = prefix_align(clip, tgt)
                    if r["t"] >= 10 and r["identity"] >= min_ident:
                        ext = cols_to_cigar(r["cols"][::-1])
                        core = ext + [(n, op) for n, op in a.cigar if op not in "SH"]
                        a = _rebuild(a, a.pos - r["t"], _merge(core), a.lclip_soft - r["q"], a.rclip_soft)
            out.append(a)
        return out

    # ................................................................ graph
    def _graph(self, orig):
        g = AssemblyGraph()
        src_len = {n: len(s) for n, s in self.project.seqs.items()}
        end_map = {}
        for w, placed in self.comps.items():
            g.segments[w] = Segment(w, self.seqs[w], len(self.seqs[w]))
            if w in self.circular or not placed:
                continue
            (o1, c1), (o2, c2) = placed[0], placed[-1]
            # working L end <- first component's outward-facing input end
            if c1[4] == "+" and c1[2] == 0:
                end_map[(c1[1], "L")] = (w, "L")
            if c1[4] == "-" and c1[3] == src_len.get(c1[1]):
                end_map[(c1[1], "R")] = (w, "L")
            if c2[4] == "+" and c2[3] == src_len.get(c2[1]):
                end_map[(c2[1], "R")] = (w, "R")
            if c2[4] == "-" and c2[2] == 0:
                end_map[(c2[1], "L")] = (w, "R")
        seen = set()
        for l in orig.links:
            ea, eb = l.ends
            wa, wb = end_map.get(ea), end_map.get(eb)
            if not wa or not wb:
                continue
            key = tuple(sorted([wa, wb]))
            if key in seen:
                continue
            seen.add(key)
            a_or = "+" if wa[1] == "R" else "-"
            b_or = "+" if wb[1] == "L" else "-"
            g.links.append(Link(wa[0], a_or, wb[0], b_or, l.overlap))
        self.end_map = end_map
        return g

    # .......................................................... hits, Hi-C
    def _map_hit(self, h):
        m0 = self.mapper.map(h.ref, h.rs)
        m1 = self.mapper.map(h.ref, h.re - 1)
        if not m0 or not m1 or m0[0] != m1[0]:
            return None
        w = m0[0]
        rs, re_ = sorted((m0[1], m1[1]))
        strand = h.strand if m0[2] == "+" else ("-" if h.strand == "+" else "+")
        return Hit(h.qname, h.qlen, h.qs, h.qe, strand, w, self.lengths[w], rs, re_ + 1, h.mapq, h.platform)

    def _map_pair(self, p):
        c1, p1, c2, p2 = p
        m1, m2 = self.mapper.map(c1, p1), self.mapper.map(c2, p2)
        if not m1 or not m2:
            return None
        return m1[0], m1[1], m2[0], m2[1]

    # ................................................ working end -> project
    def to_project_end(self, end):
        """Working contig end -> (piece, side) in the project."""
        w, side = end
        s = self.project._scaffold(self.scaffold_of[w])
        if side == "R":
            p = s["parts"][-1]
            return (p["piece"], "R" if p["orient"] == "+" else "L")
        p = s["parts"][0]
        return (p["piece"], "L" if p["orient"] == "+" else "R")

    def locate(self, w, x):
        """Working position -> (component index, piece, piece-relative cut)."""
        for i, (off, c) in enumerate(self.comps[w]):
            n = c[3] - c[2]
            if off <= x < off + n:
                piece = self.project.pieces[c[5]]
                src_cut = c[2] + (x - off) if c[4] == "+" else c[3] - (x - off)
                return i, c[5], src_cut - piece["start"], off
        return None


def _aligned(a):
    return sum(n for n, op in a.cigar if op in "M=X")


def _dedupe(alns, min_frac=0.5):
    """One read aligned to both sides of a join (e.g. primary on one old
    contig, supplementary on the other, both now extended across the
    junction) appears twice; keep the most complete copy."""
    by = {}
    for a in alns:
        by.setdefault((a.qname, a.flag & 192), []).append(a)
    out = []
    for k, lst in by.items():
        if len(lst) == 1:
            out.append(lst[0])
            continue
        lst.sort(key=lambda a: -_aligned(a))
        kept = []
        for a in lst:
            dup = False
            for b in kept:
                ov = min(a.end, b.end) - max(a.pos, b.pos)
                if a.reverse == b.reverse and ov > min_frac * (a.end - a.pos):
                    dup = True
                    break
            if not dup:
                kept.append(a)
        # the most complete copy is the primary one
        for i, a in enumerate(kept):
            sup = i > 0
            if a.supplementary != sup:
                a = replace(a, supplementary=sup, flag=(a.flag | 2048) if sup else (a.flag & ~2048))
            out.append(a)
    return out


def _merge(core):
    out = []
    for n, op in core:
        if out and out[-1][1] == op:
            out[-1] = (out[-1][0] + n, op)
        else:
            out.append((n, op))
    return out
