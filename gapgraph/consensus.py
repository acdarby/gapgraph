"""Per-column consensus with confidence, and Gap5-style problem searches.

Consensus model (in the spirit of Gap5's quality-weighted consensus): each
read base i with error probability e_i contributes weight
    w_i = ln((1 - e_i) * 3 / e_i)
to its symbol (A, C, G, T or '*' for a deletion).  The posterior of symbol b
in a column is softmax(W)[b] with W[b] = sum of weights of reads showing b,
and the consensus confidence is  Q = -10 log10(1 - P(best)), capped at 100.
e_i comes from the base quality when the alignment has one, otherwise from
a per-platform default (ONT reads usually carry qualities; simulated data
does not).

Problem searches (Gap5 contig editor "Search by ..."):
    consq     consensus quality below a threshold
    contig    contig base disagrees with the read consensus (incl. indels)
    discrep   more disagreeing reads than base qualities predict
    het       a strong second allele (mixed population / collapsed repeat copy)
    lowcov    depth below a threshold
    highcov   depth above a multiple of the median (collapsed repeats: rRNA, IS)
"""
import math
import numpy as np

SYM = "ACGT*"
_CODE = np.full(256, 255, dtype=np.uint8)
for _i, _c in enumerate("ACGT"):
    _CODE[ord(_c)] = _i
    _CODE[ord(_c.lower())] = _i

PLATFORM_ERROR = {"illumina": 0.002, "hifi": 0.002, "clr": 0.12, "ont": 0.05, "": 0.05}
LN10 = math.log(10)


def _weight(e):
    e = np.clip(e, 1e-5, 0.75)
    return np.log((1 - e) * 3 / e)


class Pileup:
    """Summary arrays for one contig."""

    def __init__(self, contig, length):
        self.contig = contig
        self.L = length
        self.W = np.zeros((length, 5), dtype=np.float64)     # weighted evidence
        self.N = np.zeros((length, 5), dtype=np.int32)       # raw counts
        self.E = np.zeros(length, dtype=np.float32)          # expected errors
        self.ins_w = np.zeros(length + 1, dtype=np.float32)  # insertion before column
        self.ins_n = np.zeros(length + 1, dtype=np.int32)

    # ------------------------------------------------------------- building
    def add_alignments(self, alns, chunk=4_000_000):
        pos_m, seg_m, len_m, err_m, qual_m = [], [], [], [], []
        pos_d, len_d, err_d = [], [], []
        pos_i, err_i = [], []
        pending = 0
        for a in alns:
            if not a.seq or a.mapq == 0 and a.supplementary:
                continue
            e_def = PLATFORM_ERROR.get(a.platform, 0.05)
            q, r = a.lclip_soft, a.pos
            for n, op in a.cigar:
                if op in "M=X":
                    pos_m.append(r)
                    seg_m.append(a.seq[q:q + n])
                    len_m.append(n)
                    if a.qual:
                        qual_m.append(a.qual[q:q + n])
                        err_m.append(-1.0)
                    else:
                        qual_m.append("")
                        err_m.append(e_def)
                    q += n
                    r += n
                    pending += n
                elif op == "I":
                    pos_i.append(r)
                    err_i.append(e_def)
                    q += n
                elif op in "DN":
                    if op == "D":
                        pos_d.append(r)
                        len_d.append(n)
                        err_d.append(e_def)
                    r += n
                    pending += n
            if pending >= chunk:
                self._flush(pos_m, seg_m, len_m, err_m, qual_m, pos_d, len_d, err_d, pos_i, err_i)
                pos_m, seg_m, len_m, err_m, qual_m = [], [], [], [], []
                pos_d, len_d, err_d = [], [], []
                pos_i, err_i = [], []
                pending = 0
        self._flush(pos_m, seg_m, len_m, err_m, qual_m, pos_d, len_d, err_d, pos_i, err_i)
        self._summarise()
        return self

    @staticmethod
    def _expand(starts, lens):
        lens = np.asarray(lens, dtype=np.int64)
        tot = int(lens.sum())
        off = np.cumsum(lens) - lens
        return np.repeat(np.asarray(starts, dtype=np.int64), lens) + (np.arange(tot) - np.repeat(off, lens))

    def _flush(self, pos_m, seg_m, len_m, err_m, qual_m, pos_d, len_d, err_d, pos_i, err_i):
        L = self.L
        if pos_m:
            refpos = self._expand(pos_m, len_m)
            codes = _CODE[np.frombuffer("".join(seg_m).encode(), dtype=np.uint8)]
            lens = np.asarray(len_m)
            e = np.repeat(np.asarray(err_m, dtype=np.float64), lens)
            if any(qual_m):
                qs = np.frombuffer("".join(x if x else "!" * n for x, n in zip(qual_m, len_m)).encode(),
                                   dtype=np.uint8).astype(np.float64) - 33
                from_q = e < 0
                e[from_q] = 10 ** (-np.maximum(qs[from_q], 1) / 10)
            ok = (codes < 4) & (refpos >= 0) & (refpos < L)
            refpos, codes, e = refpos[ok], codes[ok].astype(np.int64), e[ok]
            idx = refpos * 5 + codes
            self.W += np.bincount(idx, weights=_weight(e), minlength=L * 5).reshape(L, 5)
            self.N += np.bincount(idx, minlength=L * 5).reshape(L, 5).astype(np.int32)
            self.E += np.bincount(refpos, weights=e, minlength=L).astype(np.float32)
        if pos_d:
            refpos = self._expand(pos_d, len_d)
            e = np.repeat(np.asarray(err_d), np.asarray(len_d))
            ok = (refpos >= 0) & (refpos < L)
            refpos, e = refpos[ok], e[ok]
            self.W[:, 4] += np.bincount(refpos, weights=_weight(e), minlength=L)
            self.N[:, 4] += np.bincount(refpos, minlength=L).astype(np.int32)
            self.E += np.bincount(refpos, weights=e, minlength=L).astype(np.float32)
        if pos_i:
            p = np.clip(np.asarray(pos_i), 0, L)
            self.ins_w += np.bincount(p, weights=_weight(np.asarray(err_i)), minlength=L + 1).astype(np.float32)
            self.ins_n += np.bincount(p, minlength=L + 1).astype(np.int32)

    def _summarise(self):
        W = self.W
        self.depth = self.N.sum(axis=1)
        best = W.argmax(axis=1)
        self.cons = best.astype(np.uint8)
        mx = W.max(axis=1, keepdims=True)
        # 1 - P(best) = sum_{o != best} exp(W_o - W_best) / sum_all exp(W - W_best)
        ex = np.exp(np.clip(W - mx, -700, 0))
        tot = ex.sum(axis=1)
        p_err = np.clip((tot - 1.0) / tot, 1e-10, 1.0)
        q = -10 * np.log10(p_err)
        q[self.depth == 0] = 0
        self.qual = np.clip(q, 0, 100).astype(np.uint8)
        srt = np.sort(self.N[:, :4], axis=1)
        self.n_second = srt[:, -2]
        self.n_best_base = srt[:, -1]
        self.second = np.argsort(self.N[:, :4], axis=1)[:, -2].astype(np.uint8)
        # an insertion is in the consensus when its evidence beats the
        # evidence of the reads spanning that boundary without one
        span_w = W.sum(axis=1)
        span = np.concatenate([span_w[:1], np.minimum(span_w[:-1], span_w[1:]), span_w[-1:]])
        self.ins_in_cons = self.ins_w > 0.5 * np.maximum(span, 1e-9)
        del self.W  # summarised; keep memory per contig small

    # -------------------------------------------------------------- queries
    def consensus_string(self, start, end):
        start, end = max(0, start), min(self.L, end)
        return "".join(SYM[c] if d else " " for c, d in zip(self.cons[start:end], self.depth[start:end]))

    def masks(self, contig_seq, kind, thr):
        """Boolean array of problem columns for one search type."""
        d = self.depth
        if kind == "consq":
            return (d > 0) & (self.qual < thr)
        if kind == "contig":
            ref = _CODE[np.frombuffer(contig_seq.encode(), dtype=np.uint8)]
            ref = np.where(ref == 255, 9, ref)
            covered = d >= max(3, thr)
            return covered & ((self.cons != ref) | self.ins_in_cons[:-1])
        if kind == "discrep":
            obs = d - self.N.max(axis=1)
            return (obs >= 3) & (obs > thr * (self.E + 1.0))
        if kind == "het":
            frac = self.n_second / np.maximum(1, self.n_best_base + self.n_second)
            return (self.n_second >= 3) & (frac >= thr)
        if kind == "lowcov":
            return d < thr
        if kind == "highcov":
            med = np.median(d[d > 0]) if (d > 0).any() else 0
            return d > thr * med
        raise ValueError(f"unknown problem type {kind}")

    @staticmethod
    def regions(mask, merge=5):
        """Runs of True as [start, end), merging runs closer than `merge`."""
        if not mask.any():
            return []
        m = mask.astype(np.int8)
        edges = np.flatnonzero(np.diff(np.concatenate([[0], m, [0]])))
        runs = list(zip(edges[0::2].tolist(), edges[1::2].tolist()))
        out = [list(runs[0])]
        for s, e in runs[1:]:
            if s - out[-1][1] <= merge:
                out[-1][1] = e
            else:
                out.append([s, e])
        return out


DEFAULT_THRESHOLDS = {"consq": 30, "contig": 3, "discrep": 4.0, "het": 0.3, "lowcov": 5, "highcov": 1.8}
PROBLEM_LABELS = {
    "consq": "low consensus quality",
    "contig": "contig ≠ read consensus",
    "discrep": "excess discrepancies",
    "het": "second allele (het)",
    "lowcov": "low coverage",
    "highcov": "high coverage (collapsed repeat?)",
}
