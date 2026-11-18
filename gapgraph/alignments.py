"""Read alignments: SAM (pure Python), BAM (via optional pysam) and PAF.

Works the same for Illumina, ONT and PacBio CLR/HiFi -- the platform only
changes thresholds.  Keep soft-clips in the SAM: minimap2 `-Y` (soft-clip
supplementary alignments too), bwa mem `-Y`.
"""
import re
from bisect import bisect_left
from dataclasses import dataclass
from .seqs import open_text

_CIG = re.compile(r"(\d+)([MIDNSHP=X])")

PLATFORM_DEFAULTS = {
    #           min clip, k-mer, end window, min mapq
    "illumina": dict(min_clip=20, k=17, end_window=50, min_mapq=20),
    "hifi":     dict(min_clip=200, k=21, end_window=300, min_mapq=20),
    "clr":      dict(min_clip=500, k=13, end_window=1000, min_mapq=10),
    "ont":      dict(min_clip=500, k=13, end_window=1000, min_mapq=10),
}


@dataclass
class Alignment:
    qname: str
    ref: str
    pos: int            # 0-based leftmost aligned ref base
    end: int            # exclusive
    reverse: bool
    mapq: int
    cigar: list         # [(len, op)]
    seq: str            # SAM SEQ (ref-forward), "" if absent
    lclip: int          # soft+hard clip on the left (ref orientation)
    rclip: int
    lclip_soft: int
    rclip_soft: int
    supplementary: bool
    qlen: int           # full read length incl. hard clips
    platform: str = ""

    @property
    def lclip_seq(self):
        return self.seq[:self.lclip_soft] if self.seq else ""

    @property
    def rclip_seq(self):
        return self.seq[len(self.seq) - self.rclip_soft:] if self.seq and self.rclip_soft else ""

    def query_interval(self):
        """Aligned interval in original read coordinates (as PAF qs/qe)."""
        qs, qe = self.lclip, self.qlen - self.rclip
        if self.reverse:
            qs, qe = self.rclip, self.qlen - self.lclip
        return qs, qe


def parse_cigar(s):
    return [(int(n), op) for n, op in _CIG.findall(s)]


def _from_fields(qname, flag, ref, pos0, mapq, cigar, seq, platform):
    lclip = rclip = lsoft = rsoft = 0
    i = 0
    while i < len(cigar) and cigar[i][1] in "SH":
        lclip += cigar[i][0]
        lsoft += cigar[i][0] if cigar[i][1] == "S" else 0
        i += 1
    j = len(cigar) - 1
    while j >= i and cigar[j][1] in "SH":
        rclip += cigar[j][0]
        rsoft += cigar[j][0] if cigar[j][1] == "S" else 0
        j -= 1
    rlen = sum(n for n, op in cigar if op in "MDN=X")
    qlen = sum(n for n, op in cigar if op in "MIS=XH")
    return Alignment(qname, ref, pos0, pos0 + rlen, bool(flag & 16), mapq, cigar,
                     "" if seq == "*" else seq, lclip, rclip, lsoft, rsoft,
                     bool(flag & 2048), qlen, platform)


def read_sam(path, platform=""):
    refs = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    t = dict(x.split(":", 1) for x in line.rstrip().split("\t")[1:])
                    refs[t["SN"]] = int(t["LN"])
                continue
            f = line.rstrip("\n").split("\t")
            flag = int(f[1])
            if flag & 4 or flag & 256 or f[5] == "*":
                continue
            yield _from_fields(f[0], flag, f[2], int(f[3]) - 1, int(f[4]), parse_cigar(f[5]), f[9], platform)


_PYSAM_OPS = "MIDNSHP=X"


def read_bam(path, platform="", region=None):
    import pysam  # optional dependency
    with pysam.AlignmentFile(path) as bam:
        it = bam.fetch(*region) if region else bam.fetch(until_eof=True)
        for r in it:
            if r.is_unmapped or r.is_secondary:
                continue
            cig = [(n, _PYSAM_OPS[op]) for op, n in r.cigartuples]
            yield _from_fields(r.query_name, r.flag, r.reference_name, r.reference_start,
                               r.mapping_quality, cig, r.query_sequence or "*", platform)


def read_alignments(path, platform=""):
    p = str(path)
    if p.endswith(".bam") or p.endswith(".cram"):
        return read_bam(p, platform)
    return read_sam(p, platform)


class AlignmentStore:
    """In-memory per-contig index (fine for bacterial / BAC / region-sized work;
    for whole large genomes use indexed BAM + pysam region fetches)."""

    def __init__(self):
        self.by_ref = {}
        self._starts = {}
        self._maxlen = {}

    def add(self, alns):
        for a in alns:
            self.by_ref.setdefault(a.ref, []).append(a)
        for r, lst in self.by_ref.items():
            lst.sort(key=lambda a: a.pos)
            self._starts[r] = [a.pos for a in lst]
            self._maxlen[r] = max((a.end - a.pos for a in lst), default=0)
        return self

    def all(self):
        for lst in self.by_ref.values():
            yield from lst

    def fetch(self, ref, start, end):
        lst = self.by_ref.get(ref, [])
        lo = bisect_left(self._starts.get(ref, []), start - self._maxlen.get(ref, 0))
        for a in lst[lo:]:
            if a.pos >= end:
                break
            if a.end + a.rclip_soft > start and a.pos - a.lclip_soft < end:
                yield a

    def coverage(self, ref, length, bins=200):
        step = max(1, length // bins)
        cov = [0] * (length // step + 1)
        for a in self.by_ref.get(ref, []):
            for b in range(a.pos // step, min(len(cov), a.end // step + 1)):
                cov[b] += 1
        return step, cov


@dataclass
class Hit:
    """Read-to-contig hit in PAF terms (from PAF or derived from SAM)."""
    qname: str
    qlen: int
    qs: int
    qe: int
    strand: str
    ref: str
    rlen: int
    rs: int
    re: int
    mapq: int
    platform: str = ""


def read_paf(path, platform=""):
    with open_text(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 12:
                continue
            if any(t == "tp:A:S" for t in f[12:]):
                continue  # secondary
            yield Hit(f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]),
                      int(f[7]), int(f[8]), int(f[11]), platform)


def hits_from_alignments(alns, ref_lengths):
    for a in alns:
        qs, qe = a.query_interval()
        yield Hit(a.qname, a.qlen, qs, qe, "-" if a.reverse else "+", a.ref,
                  ref_lengths.get(a.ref, a.end), a.pos, a.end, a.mapq, a.platform)


def render_read(a, start, end):
    """Render an alignment in reference coordinates for the editor.

    Returns dict with `x` (ref coord of first char), `s` (one char per ref
    column: bases upper-case, soft-clipped bases lower-case, '*' deletion) and
    `ins` (ref columns that have an insertion before them, with the bases).
    """
    chars, ins = [], []
    x0 = a.pos - a.lclip_soft
    q, r = 0, a.pos
    if a.lclip_soft:
        chars.extend(a.seq[:a.lclip_soft].lower() if a.seq else "?" * a.lclip_soft)
        q = a.lclip_soft
    for n, op in a.cigar:
        if op in "M=X":
            chars.extend(a.seq[q:q + n] if a.seq else "=" * n)
            q += n
            r += n
        elif op == "I":
            ins.append([r, a.seq[q:q + n] if a.seq else str(n)])
            q += n
        elif op in "DN":
            chars.extend("*" * n)
            r += n
    if a.rclip_soft:
        chars.extend(a.seq[q:q + a.rclip_soft].lower() if a.seq else "?" * a.rclip_soft)
    s = "".join(chars)
    lo, hi = max(start, x0), min(end, x0 + len(s))
    return {"x": lo, "s": s[lo - x0:hi - x0], "ins": [i for i in ins if lo <= i[0] < hi],
            "rev": a.reverse, "name": a.qname, "mapq": a.mapq, "supp": a.supplementary,
            "platform": a.platform, "pos": a.pos, "end": a.end}
