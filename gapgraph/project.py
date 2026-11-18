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
        self._undo.append(copy.deepcopy((self.pieces, self.scaffolds, self.discarded, self._counter)))
        self.log.append({"t": time.strftime("%Y-%m-%d %H:%M:%S"), "op": op, **kw})

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
        if name in self.pieces and not self.pieces[name].get("copy"):
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
        parts = s["parts"][::-1]
        gaps = [p["gap"] for p in s["parts"]][:-1][::-1] + [0]
        s["parts"] = [{"piece": p["piece"], "orient": flip_orient(p["orient"]), "gap": g}
                      for p, g in zip(parts, gaps)]

    def join(self, end_a, end_b, gap=100, reason="", evidence=None, via=()):
        sa, side_a = self._scaffold_end(tuple(end_a))
        sb, side_b = self._scaffold_end(tuple(end_b))
        if sa == sb:
            raise ProjectError("both ends on the same scaffold (circularisation not supported here)")
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

    # ----------------------------------------------------------------- export
    def scaffold_seq(self, s):
        out = []
        trim = 0
        for part in s["parts"]:
            seq = self.piece_seq(part["piece"])
            if part["orient"] == "-":
                seq = revcomp(seq)
            out.append(seq[trim:])
            g = part["gap"]
            if g > 0:
                out.append("N" * g)
                trim = 0
            else:
                trim = -g  # overlap: trim from the next piece
        return "".join(out)

    def export(self, prefix):
        recs = [(s["name"], self.scaffold_seq(s)) for s in self.scaffolds]
        write_fasta(recs, prefix + ".fasta")
        with open(prefix + ".agp", "w") as fh:
            fh.write("##agp-version\t2.1\n")
            for s in self.scaffolds:
                pos, k, trim = 1, 1, 0
                for i, part in enumerate(s["parts"]):
                    p = self.pieces[part["piece"]]
                    ln = p["end"] - p["start"] - trim
                    cs, ce = (p["start"] + 1 + trim, p["end"]) if part["orient"] == "+" else (p["start"] + 1, p["end"] - trim)
                    fh.write(f"{s['name']}\t{pos}\t{pos + ln - 1}\t{k}\tW\t{p['src']}\t{cs}\t{ce}\t{part['orient']}\n")
                    pos += ln
                    k += 1
                    g = part["gap"]
                    trim = -g if g < 0 else 0
                    if g > 0 and i < len(s["parts"]) - 1:
                        fh.write(f"{s['name']}\t{pos}\t{pos + g - 1}\t{k}\tN\t{g}\tscaffold\tyes\t{part.get('ev', 'unspecified')}\n")
                        pos += g
                        k += 1
        with open(prefix + ".log.json", "w") as fh:
            json.dump(self.log, fh, indent=1)
        return [prefix + ".fasta", prefix + ".agp", prefix + ".log.json"]
