"""Minimal GFA1 assembly graph with contig *ends* as first-class objects.

A contig end is (name, "L") or (name, "R").  A GFA link a(+/-) -> b(+/-) joins:
    a+ -> b+ : a.R -- b.L        a+ -> b- : a.R -- b.R
    a- -> b+ : a.L -- b.L        a- -> b- : a.L -- b.R
Every join in gapgraph (graph link, soft-clip, long read, Hi-C) is an unordered
pair of ends, so evidence from different sources lands on the same key.
"""
from dataclasses import dataclass, field
from .seqs import open_text


def join_key(end_a, end_b):
    return tuple(sorted([tuple(end_a), tuple(end_b)]))


def end_from_orient(name, orient, outgoing):
    """End used when leaving (outgoing=True) or entering an oriented segment."""
    if outgoing:
        return (name, "R" if orient == "+" else "L")
    return (name, "L" if orient == "+" else "R")


def other_end(end):
    return (end[0], "L" if end[1] == "R" else "R")


@dataclass
class Segment:
    name: str
    seq: str
    length: int
    tags: dict = field(default_factory=dict)


@dataclass
class Link:
    a: str
    a_or: str
    b: str
    b_or: str
    overlap: int = 0

    @property
    def ends(self):
        return end_from_orient(self.a, self.a_or, True), end_from_orient(self.b, self.b_or, False)


def _cigar_overlap(cig):
    if cig in ("*", "", None):
        return 0
    n, num = 0, ""
    for ch in cig:
        if ch.isdigit():
            num += ch
        else:
            if ch in "MX=D":
                n += int(num or 0)
            num = ""
    return n


class AssemblyGraph:
    def __init__(self):
        self.segments = {}
        self.links = []
        self.paths = {}

    # ---------------------------------------------------------------- io
    @classmethod
    def from_gfa(cls, path):
        g = cls()
        with open_text(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == "S":
                    seq = f[2] if f[2] != "*" else ""
                    tags = {t[:2]: t[5:] for t in f[3:] if len(t) > 5}
                    length = len(seq) or int(tags.get("LN", 0))
                    g.segments[f[1]] = Segment(f[1], seq, length, tags)
                elif f[0] == "L":
                    g.links.append(Link(f[1], f[2], f[3], f[4], _cigar_overlap(f[5] if len(f) > 5 else "*")))
                elif f[0] == "P":
                    g.paths[f[1]] = [(s[:-1], s[-1]) for s in f[2].split(",")]
        return g

    @classmethod
    def from_fasta(cls, seqs):
        g = cls()
        for n, s in seqs.items():
            g.segments[n] = Segment(n, s, len(s))
        return g

    def to_gfa(self, path, extra_paths=None):
        with open(path, "w") as fh:
            fh.write("H\tVN:Z:1.0\n")
            for s in self.segments.values():
                fh.write(f"S\t{s.name}\t{s.seq or '*'}\tLN:i:{s.length}\n")
            for l in self.links:
                fh.write(f"L\t{l.a}\t{l.a_or}\t{l.b}\t{l.b_or}\t{l.overlap}M\n")
            for name, p in {**self.paths, **(extra_paths or {})}.items():
                fh.write(f"P\t{name}\t{','.join(n + o for n, o in p)}\t*\n")

    # ------------------------------------------------------------ queries
    def link_keys(self):
        return {join_key(*l.ends): l for l in self.links}

    def end_degree(self):
        deg = {}
        for l in self.links:
            for e in l.ends:
                deg[e] = deg.get(e, 0) + 1
        return deg

    def dead_ends(self):
        deg = self.end_degree()
        return [(n, side) for n in self.segments for side in "LR" if deg.get((n, side), 0) == 0]

    def tangles(self):
        """Segments with >1 link on both ends: collapsed repeats needing resolution."""
        deg = self.end_degree()
        return [n for n in self.segments if deg.get((n, "L"), 0) > 1 and deg.get((n, "R"), 0) > 1]

    def to_json(self):
        return {
            "segments": [{"name": s.name, "length": s.length} for s in self.segments.values()],
            "links": [{"a": list(l.ends[0]), "b": list(l.ends[1]), "overlap": l.overlap} for l in self.links],
            "dead_ends": [list(e) for e in self.dead_ends()],
            "tangles": self.tangles(),
        }
