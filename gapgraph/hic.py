"""Hi-C contacts summarised at contig-end resolution.

Contacts decay with distance, so for two contigs that are neighbours the
end-pair that is physically adjacent collects most contacts: this gives
orientation (which ends join) as well as adjacency.  Input: 4DN `.pairs`
(readID chr1 pos1 chr2 pos2 ...) e.g. from pairtools / chromap / juicer
(convert merged_nodups with a one-liner), on contig coordinates.
"""
import math
from collections import defaultdict
from .seqs import open_text
from .gfa import join_key


def read_pairs(path):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 5:
                yield f[1], int(f[2]) - 1, f[3], int(f[4]) - 1


class HiC:
    def __init__(self, contig_len, end_bin=20000, bins=40):
        self.L = contig_len
        self.end_bin = end_bin
        self.end_counts = defaultdict(int)     # join_key -> count
        self.contig_counts = defaultdict(int)  # (c1,c2) sorted -> count
        self.totals = defaultdict(int)         # contig -> contacts with other contigs
        self.intra = defaultdict(int)
        self.end_totals = defaultdict(int)     # end -> inter-contig end-end contacts
        self.n = 0

    def _end(self, c, p):
        L = self.L[c]
        w = min(self.end_bin, L // 2)
        if p < w:
            return (c, "L")
        if p >= L - w:
            return (c, "R")
        return None

    def add(self, pairs):
        for c1, p1, c2, p2 in pairs:
            if c1 not in self.L or c2 not in self.L:
                continue
            self.n += 1
            if c1 == c2:
                self.intra[c1] += 1
                continue
            self.contig_counts[tuple(sorted((c1, c2)))] += 1
            self.totals[c1] += 1
            self.totals[c2] += 1
            e1, e2 = self._end(c1, p1), self._end(c2, p2)
            if e1 and e2:
                self.end_counts[join_key(e1, e2)] += 1
                self.end_totals[e1] += 1
                self.end_totals[e2] += 1
        return self

    def score(self, key):
        (a, sa), (b, sb) = key
        if a == b:
            return None
        n = self.end_counts.get(key, 0)
        four = sum(self.end_counts.get(join_key((a, x), (b, y)), 0) for x in "LR" for y in "LR")
        # share: fraction of each end's inter-contig contacts that go to the other end
        share = math.sqrt(n / max(1, self.end_totals[(a, sa)]) * n / max(1, self.end_totals[(b, sb)]))
        return {"contacts": n, "norm": round(share, 4),
                "orient_frac": round(n / four, 3) if four else 0.0}

    def candidate_keys(self, top_per_end=3, min_contacts=5):
        best = defaultdict(list)
        for key, n in self.end_counts.items():
            if n >= min_contacts:
                for e in key:
                    best[e].append((n, key))
        keys = set()
        for e, lst in best.items():
            for n, k in sorted(lst, reverse=True)[:top_per_end]:
                keys.add(k)
        return keys

    def matrix(self, order):
        """log10 normalised contig-level matrix for the heatmap."""
        m = []
        for a in order:
            row = []
            for b in order:
                if a == b:
                    v = self.intra.get(a, 0) / max(1, self.L[a]) * 1000
                else:
                    v = self.contig_counts.get(tuple(sorted((a, b))), 0) / math.sqrt(
                        max(1, self.L[a]) * max(1, self.L[b])) * 1000
                row.append(round(math.log10(1 + v), 3))
            m.append(row)
        return {"order": order, "m": m}
