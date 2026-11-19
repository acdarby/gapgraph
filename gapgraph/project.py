"""Editable finishing project: scaffolds over (sub)contigs, with undo and an
audit log -- the Gap5 'database' idea, kept as a JSON file next to the inputs.
Alignments stay in their (indexed) files; the project only stores layout."""
import copy
import json
import time
from .seqs import revcomp, write_fasta


class ProjectError(Exception):
    pass


def flip_orient(o):
    return "-" if o == "+" else "+"


class Project:
    def __init__(self, contig_seqs, path=None):
        self.seqs = contig_seqs
        self.path = path
        self.pieces = {n: {"src": n, "start": 0, "end": len(s)} for n, s in contig_seqs.items()}
        self.scaffolds = [{"name": f"scf_{n}", "parts": [{"piece": n, "orient": "+", "gap": 0}]}
                          for n in contig_seqs]
        self.log = []
        self.discarded = []
        self._undo = []
        self._counter = 0

    # ------------------------------------------------------------- persistence
    def state(self):
        return {"pieces": self.pieces, "scaffolds": self.scaffolds, "discarded": self.discarded,
                "log": self.log, "counter": self._counter}

    def save(self):
        if self.path:
            with open(self.path, "w") as fh:
                json.dump(self.state(), fh, indent=1)

    def load(self):
        with open(self.path) as fh:
            st = json.load(fh)
        self.pieces, self.scaffolds, self.log = st["pieces"], st["scaffolds"], st["log"]
        self._counter = st.get("counter", 0)
        self.discarded = st.get("discarded", [])

    def _checkpoint(self, op, **kw):
        if not getattr(self, "_tx", 0):
            self._undo.append(copy.deepcopy((self.pieces, self.scaffolds, self.discarded, self._counter)))
        self.log.append({"t": time.strftime("%Y-%m-%d %H:%M:%S"), "op": op, **kw})

    def transaction(self, op, **kw):
        """Group several operations so that one undo reverses them all."""
        proj = self

        class _Tx:
            def __enter__(self_):
                proj._checkpoint(op, **kw)
                proj._tx = getattr(proj, "_tx", 0) + 1

            def __exit__(self_, *exc):
                proj._tx -= 1
                if exc[0] is not None and not proj._tx:   # roll back a failed transaction
                    proj.pieces, proj.scaffolds, proj.discarded, proj._counter = proj._undo.pop()
                proj.save()
                return False
        return _Tx()

    def undo(self):
        if not self._undo:
            raise ProjectError("nothing to undo")
        self.pieces, self.scaffolds, self.discarded, self._counter = self._undo.pop()
        self.log.append({"t": time.strftime("%Y-%m-%d %H:%M:%S"), "op": "undo"})
        self.save()

    # ---------------------------------------------------------------- lookups
    def piece_seq(self, name):
        p = self.pieces[name]
        return self.seqs[p["src"]][p["start"]:p["end"]]

    def piece_for_end(self, end):
        """Map an original contig end to the piece that currently carries it."""
        c, side = end
        for n, p in self.pieces.items():
            if p.get("copy"):
                continue
            if p["src"] == c and ((side == "L" and p["start"] == 0) or
                                  (side == "R" and p["end"] == len(self.seqs[c]))):
                return n
        raise ProjectError(f"end {c}:{side} no longer exists")

    def find(self, piece):
        for si, s in enumerate(self.scaffolds):
            for pi, part in enumerate(s["parts"]):
                if part["piece"] == piece:
                    return si, pi
        raise ProjectError(f"piece {piece} not in any scaffold")

    def resolve_end(self, end):
        """[piece_or_contig, side] -> (piece, side relative to the piece)."""
        name, side = end
        if name in self.pieces:
            return name, side
        return self.piece_for_end((name, side)), side

    def _scaffold_end(self, end):
        """Return (scaffold index, 'L'|'R') if this contig end is a scaffold terminus."""
        piece, _ = self.resolve_end(end)
        si, pi = self.find(piece)
        s = self.scaffolds[si]
        part = s["parts"][pi]
        # which side of the scaffold does this end of the piece face?
        faces_right = (end[1] == "R") == (part["orient"] == "+")
        if faces_right and pi == len(s["parts"]) - 1:
            return si, "R"
        if not faces_right and pi == 0:
            return si, "L"
        raise ProjectError(f"{end[0]}:{end[1]} is internal to scaffold {s['name']}")

    # ------------------------------------------------------------- operations
    def _reverse(self, s):
        if s.get("circular") and s.get("origin"):
            # unrotated circle = linear[:U - t]; reversed it is revcomp(linear)[:U - t]
            U = len(self.scaffold_seq({**s, "circular": False}, rotate=False))
            t = max(0, -s.get("circ_gap", 0))
            o = s["origin"]
            s["origin"] = {"pos": (U - 1 - o["pos"]) % (U - t), "strand": flip_orient(o["strand"])}
        parts = s["parts"][::-1]
        gaps = [(p["gap"], p.get("ev", "unspecified")) for p in s["parts"]][:-1][::-1] + [(0, "unspecified")]
        s["parts"] = [{"piece": p["piece"], "orient": flip_orient(p["orient"]), "gap": g, "ev": ev}
                      for p, (g, ev) in zip(parts, gaps)]

    def join(self, end_a, end_b, gap=100, reason="", evidence=None, via=()):
        sa, side_a = self._scaffold_end(tuple(end_a))
        sb, side_b = self._scaffold_end(tuple(end_b))
        if sa == sb:
            if side_a == side_b:
                raise ProjectError("same end twice")
            return self.circularise(self.scaffolds[sa]["name"], gap, reason=reason)
        self._checkpoint("join", a=list(end_a), b=list(end_b), gap=gap, reason=reason)
        A, B = self.scaffolds[sa], self.scaffolds[sb]
        if side_a == "L":
            self._reverse(A)
        if side_b == "R":
            self._reverse(B)
        ev = "proximity_ligation" if evidence == ["hic"] else "unspecified"
        for v in via:  # insert copies of a collapsed repeat between the two ends
            name, orient = v[:-1], v[-1]
            self._counter += 1
            copy_name = f"{name}#{self._counter}"
            self.pieces[copy_name] = {"src": name, "start": 0, "end": len(self.seqs[name]), "copy": True}
            A["parts"][-1]["gap"] = 0
            A["parts"][-1]["ev"] = ev
            A["parts"].append({"piece": copy_name, "orient": orient, "gap": 0})
        A["parts"][-1]["gap"] = int(gap)
        A["parts"][-1]["ev"] = ev
        A["parts"].extend(B["parts"])
        A["name"] = f"{A['name']}+{B['name'].removeprefix('scf_')}"[:80]
        self.scaffolds.pop(sb)
        self.save()
        return A["name"]

    def break_scaffold(self, scaffold, after_index):
        si = next(i for i, s in enumerate(self.scaffolds) if s["name"] == scaffold)
        s = self.scaffolds[si]
        if not 0 <= after_index < len(s["parts"]) - 1:
            raise ProjectError("break index out of range")
        self._checkpoint("break", scaffold=scaffold, after=after_index)
        right = s["parts"][after_index + 1:]
        s["parts"] = s["parts"][:after_index + 1]
        s["parts"][-1]["gap"] = 0
        self._counter += 1
        self.scaffolds.insert(si + 1, {"name": f"{scaffold}.b{self._counter}", "parts": right})
        self.save()

    def split_contig(self, piece, pos, reason=""):
        """Break a contig at a (piece-relative) position, e.g. a clip hotspot."""
        p = self.pieces[piece]
        if not 0 < pos < p["end"] - p["start"]:
            raise ProjectError("split position outside contig")
        self._checkpoint("split", piece=piece, pos=pos, reason=reason)
        cut = p["start"] + pos
        left = f"{p['src']}:{p['start']}-{cut}"
        right = f"{p['src']}:{cut}-{p['end']}"
        self.pieces.pop(piece)
        self.pieces[left] = {"src": p["src"], "start": p["start"], "end": cut}
        self.pieces[right] = {"src": p["src"], "start": cut, "end": p["end"]}
        si, pi = self.find(piece)
        s = self.scaffolds[si]
        part = s["parts"][pi]
        a, b = (left, right) if part["orient"] == "+" else (right, left)
        before, after = s["parts"][:pi], s["parts"][pi + 1:]
        s["parts"] = before + [{"piece": a, "orient": part["orient"], "gap": 0}]
        self._counter += 1
        new = {"name": f"{s['name']}.s{self._counter}",
               "parts": [{"piece": b, "orient": part["orient"], "gap": part["gap"]}] + after}
        self.scaffolds.insert(si + 1, new)
        self.save()

    def discard(self, scaffold, reason=""):
        """Move a scaffold (e.g. a chimeric fragment) out of the assembly."""
        si = next(i for i, s in enumerate(self.scaffolds) if s["name"] == scaffold)
        self._checkpoint("discard", scaffold=scaffold, reason=reason)
        self.discarded.append(self.scaffolds.pop(si))
        self.save()

    def flip(self, scaffold):
        s = next(s for s in self.scaffolds if s["name"] == scaffold)
        self._checkpoint("flip", scaffold=scaffold)
        self._reverse(s)
        self.save()

    def circularise(self, scaffold, gap=0, reason=""):
        """Close a scaffold into a circle.  gap < 0: the two ends overlap by
        -gap bp (duplicated sequence, typical of Canu/Flye/hifiasm circular
        contigs) and that much is trimmed from the end; gap > 0: Ns inserted."""
        s = self._scaffold(scaffold)
        if s.get("circular"):
            raise ProjectError(f"{scaffold} is already circular")
        if -gap >= len(self.scaffold_seq(s, rotate=False)):
            raise ProjectError("overlap longer than the scaffold")
        self._checkpoint("circularise", scaffold=scaffold, gap=gap, reason=reason)
        s["circular"] = True
        s["circ_gap"] = int(gap)
        self.save()
        return scaffold

    def linearise(self, scaffold):
        s = self._scaffold(scaffold)
        self._checkpoint("linearise", scaffold=scaffold)
        for k in ("circular", "circ_gap", "origin"):
            s.pop(k, None)
        self.save()

    def set_origin(self, scaffold, pos, strand="+", reason=""):
        """Rotate a circular scaffold so `pos` (0-based, in the unrotated
        circular sequence) becomes base 1; strand '-' also reverse-complements
        so a minus-strand dnaA ends up forward.  `pos` is the gene's first base."""
        s = self._scaffold(scaffold)
        if not s.get("circular"):
            raise ProjectError("circularise the scaffold before setting its start")
        L = len(self.scaffold_seq(s, rotate=False))
        if not 0 <= pos < L:
            raise ProjectError("start position outside scaffold")
        self._checkpoint("set_origin", scaffold=scaffold, pos=pos, strand=strand, reason=reason)
        s["origin"] = {"pos": int(pos), "strand": strand}
        self.save()

    def _scaffold(self, name):
        for s in self.scaffolds:
            if s["name"] == name:
                return s
        raise ProjectError(f"no scaffold {name}")

    # ----------------------------------------------------------------- export
    def components(self, s, rotate=True):
        """Scaffold as AGP-like components: ('W', src, start, end, orient) with
        0-based half-open src coords, or ('N', length, evidence)."""
        comps, trim = [], 0
        for i, part in enumerate(s["parts"]):
            p = self.pieces[part["piece"]]
            a, b = p["start"], p["end"]
            if part["orient"] == "+":
                a += trim
            else:
                b -= trim
            comps.append(["W", p["src"], a, b, part["orient"], part["piece"]])
            g = part["gap"]
            trim = 0
            if i < len(s["parts"]) - 1:
                if g > 0:
                    comps.append(["N", g, part.get("ev", "unspecified")])
                else:
                    trim = -g
        if s.get("circular"):
            g = s.get("circ_gap", 0)
            if g > 0:
                comps.append(["N", g, "unspecified"])
            t = -g
            while t > 0:          # trim the duplicated overlap off the end
                c = comps[-1]
                n = c[1] if c[0] == "N" else c[3] - c[2]
                if n <= t:
                    comps.pop()
                    t -= n
                    continue
                if c[0] == "N":
                    c[1] -= t
                elif c[4] == "+":
                    c[3] -= t
                else:
                    c[2] += t
                t = 0
            if rotate and s.get("origin"):
                comps = self._rotate(comps, s["origin"])
        return comps

    @staticmethod
    def _clen(c):
        return c[1] if c[0] == "N" else c[3] - c[2]

    def _rotate(self, comps, origin):
        L = sum(self._clen(c) for c in comps)
        pos = origin["pos"] % L
        if origin.get("strand", "+") == "-":
            comps = [c if c[0] == "N" else ["W", c[1], c[2], c[3], flip_orient(c[4]), c[5]] for c in comps[::-1]]
            pos = L - 1 - pos
        out, x = [], 0
        for i, c in enumerate(comps):
            n = self._clen(c)
            if x <= pos < x + n:
                k = pos - x
                if k == 0:
                    return comps[i:] + comps[:i]
                if c[0] == "N":
                    left, right = ["N", k, c[2]], ["N", n - k, c[2]]
                elif c[4] == "+":
                    left, right = ["W", c[1], c[2], c[2] + k, "+", c[5]], ["W", c[1], c[2] + k, c[3], "+", c[5]]
                else:
                    left, right = ["W", c[1], c[3] - k, c[3], "-", c[5]], ["W", c[1], c[2], c[3] - k, "-", c[5]]
                return [right] + comps[i + 1:] + comps[:i] + [left]
            x += n
        return comps

    def scaffold_seq(self, s, rotate=True):
        out = []
        for c in self.components(s, rotate):
            if c[0] == "N":
                out.append("N" * c[1])
            else:
                seq = self.seqs[c[1]][c[2]:c[3]]
                out.append(seq if c[4] == "+" else revcomp(seq))
        return "".join(out)

    def export(self, prefix):
        recs = []
        for s in self.scaffolds:
            seq = self.scaffold_seq(s)
            name = s["name"]
            if s.get("circular"):
                name += f" length={len(seq)} circular=true [topology=circular]"
            recs.append((name, seq))
        write_fasta(recs, prefix + ".fasta")
        with open(prefix + ".agp", "w") as fh:
            fh.write("##agp-version\t2.1\n")
            for s in self.scaffolds:
                if s.get("circular"):
                    fh.write(f"# {s['name']} is circular\n")
                pos, k = 1, 1
                for c in self.components(s):
                    n = self._clen(c)
                    if c[0] == "W":
                        fh.write(f"{s['name']}\t{pos}\t{pos + n - 1}\t{k}\tW\t{c[1]}\t{c[2] + 1}\t{c[3]}\t{c[4]}\n")
                    else:
                        fh.write(f"{s['name']}\t{pos}\t{pos + n - 1}\t{k}\tN\t{n}\tscaffold\tyes\t{c[2]}\n")
                    pos += n
                    k += 1
        with open(prefix + ".log.json", "w") as fh:
            json.dump(self.log, fh, indent=1)
        return [prefix + ".fasta", prefix + ".agp", prefix + ".log.json"]
