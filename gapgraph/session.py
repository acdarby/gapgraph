"""Load inputs, compute all evidence, answer editor queries."""
import os
from .gfa import AssemblyGraph
from .seqs import read_fasta
from .alignments import (AlignmentStore, read_alignments, read_paf, hits_from_alignments,
                         render_read, PLATFORM_DEFAULTS)
from . import softclip, longreads, evidence
from .hic import HiC, read_pairs
from .project import Project


class Session:
    def __init__(self, gfa=None, fasta=None, alignments=(), paf=(), pairs=None, project=None,
                 hic_end_bin=20000):
        """alignments / paf: lists of (path, platform) with platform in
        illumina | hifi | clr | ont."""
        if gfa:
            self.graph = AssemblyGraph.from_gfa(gfa)
            if fasta:
                for n, s in read_fasta(fasta).items():
                    if n in self.graph.segments:
                        self.graph.segments[n].seq = s
        else:
            self.graph = AssemblyGraph.from_fasta(read_fasta(fasta))
        self.seqs = {n: s.seq for n, s in self.graph.segments.items()}
        missing = [n for n, s in self.seqs.items() if not s]
        if missing:
            raise ValueError(f"no sequence for segments {missing[:5]} -- pass --fasta")
        self.lengths = {n: len(s) for n, s in self.seqs.items()}

        self.stores = {}
        for path, plat in alignments:
            self.stores.setdefault(plat, AlignmentStore()).add(read_alignments(path, plat))

        self.hits = []
        for path, plat in paf:
            self.hits.extend(read_paf(path, plat))
        for plat, st in self.stores.items():
            if plat != "illumina":
                self.hits.extend(hits_from_alignments(st.all(), self.lengths))

        self.hic = None
        if pairs:
            self.hic = HiC(self.lengths, end_bin=hic_end_bin).add(read_pairs(pairs))

        self.project = Project(self.seqs, path=project)
        if project and os.path.exists(project):
            self.project.load()
        self.compute()

    def compute(self):
        clip_joins, self.hotspots = {}, []
        self.clips = {}
        for plat, st in self.stores.items():
            d = PLATFORM_DEFAULTS.get(plat, PLATFORM_DEFAULTS["ont"])
            clips = softclip.outward_clips(st, self.lengths, d["min_clip"], d["end_window"], d["min_mapq"])
            for e, lst in clips.items():
                self.clips.setdefault(e, []).extend(lst)
            for key, v in softclip.match_clips(clips, self.seqs, d["k"]).items():
                m = clip_joins.setdefault(key, {"reads": 0, "gap": v["gap"], "platforms": {}})
                m["reads"] += v["reads"]
                m["platforms"].update(v["platforms"])
            for h in softclip.internal_hotspots(st, self.lengths, d["min_clip"], d["end_window"], d["min_mapq"]):
                h["platform"] = plat
                self.hotspots.append(h)
        self.lr_joins, self.lr_paths = longreads.bridges(self.hits)
        self.candidates = evidence.build_candidates(self.graph, clip_joins, self.lr_joins, self.hic, self.lr_paths)
        self.cand_by_id = {c["id"]: c for c in self.candidates}

    # ------------------------------------------------------------ API views
    def overview(self):
        g = self.graph.to_json()
        for s in g["segments"]:
            cov = {}
            for plat, st in self.stores.items():
                n = sum(a.end - a.pos for a in st.by_ref.get(s["name"], []))
                cov[plat] = round(n / max(1, s["length"]), 1)
            s["cov"] = cov
            s["clips"] = {side: len(self.clips.get((s["name"], side), [])) for side in "LR"}
        return {
            "graph": g,
            "candidates": self.candidates,
            "hotspots": self.hotspots,
            "paths": self.lr_paths[:200],
            "scaffolds": self.project.scaffolds,
            "discarded": self.project.discarded,
            "pieces": self.project.pieces,
            "log": self.project.log[-50:],
            "platforms": sorted(self.stores),
            "has_hic": self.hic is not None,
        }

    def region(self, contig, start, end, max_rows=400, overhang=20000):
        """Editor window.  The window may extend past either contig end so the
        soft-clipped bases of reads running off the end stay visible."""
        L = self.lengths[contig]
        start, end = max(-overhang, start), min(L + overhang, end)
        rows, row_end, reads = [], [], []
        for plat in sorted(self.stores):
            reads.extend(self.stores[plat].fetch(contig, start, end))
        reads.sort(key=lambda a: a.pos - a.lclip_soft)
        for a in reads:
            r = render_read(a, start, end)
            if not r["s"]:
                continue
            x0 = r["x"]
            for i, e in enumerate(row_end):
                if e <= x0:
                    rows[i].append(r)
                    row_end[i] = x0 + len(r["s"]) + 1
                    break
            else:
                if len(rows) >= max_rows:
                    continue
                rows.append([r])
                row_end.append(x0 + len(r["s"]) + 1)
        depth = [0] * (end - start)
        for a in reads:
            for x in range(max(start, a.pos), min(end, a.end)):
                depth[x - start] += 1
        ref = " " * max(0, -start) + self.seqs[contig][max(0, start):min(L, end)] + " " * max(0, end - L)
        return {"contig": contig, "start": start, "end": end, "length": L,
                "ref": ref, "rows": rows, "depth": depth,
                "consensus": self._consensus(rows, start, end)}

    @staticmethod
    def _consensus(rows, start, end):
        from collections import Counter
        cols = [Counter() for _ in range(end - start)]
        for row in rows:
            for r in row:
                for i, ch in enumerate(r["s"]):
                    if ch.isupper() or ch == "*":
                        cols[r["x"] - start + i][ch] += 1
        return "".join(c.most_common(1)[0][0] if c else " " for c in cols)

    def end_detail(self, contig, side):
        clips = self.clips.get((contig, side), [])
        seqs = [c for _, c, _ in clips]
        return {"end": [contig, side], "n_clips": len(clips),
                "clip_lengths": sorted((len(s) for s in seqs), reverse=True)[:200],
                "consensus": softclip.naive_consensus(seqs),
                "candidates": [c for c in self.candidates if [contig, side] in (c["a"], c["b"])]}

    def hic_matrix(self):
        if not self.hic:
            return None
        order = []
        for s in self.project.scaffolds:
            for p in s["parts"]:
                src = self.project.pieces[p["piece"]]["src"]
                if src not in order:
                    order.append(src)
        return self.hic.matrix(order)

    def op(self, req):
        P = self.project
        op = req.get("op")
        if op == "join":
            c = self.cand_by_id.get(req.get("id"))
            a, b = (c["a"], c["b"]) if c else (req["a"], req["b"])
            gap = int(req.get("gap", c["gap"] if c else 100))
            return {"scaffold": P.join(a, b, gap, reason=req.get("reason", ""),
                                       evidence=c["evidence"] if c else None,
                                       via=req.get("via", c["via"] if c else []))}
        if op == "break":
            return P.break_scaffold(req["scaffold"], int(req["after"]))
        if op == "split":
            return P.split_contig(req["piece"], int(req["pos"]), reason=req.get("reason", ""))
        if op == "discard":
            return P.discard(req["scaffold"], reason=req.get("reason", ""))
        if op == "flip":
            return P.flip(req["scaffold"])
        if op == "undo":
            return P.undo()
        if op == "auto":
            return self.auto_scaffold(float(req.get("min_score", 6)))
        raise ValueError(f"unknown op {op}")

    def auto_scaffold(self, min_score=6.0):
        """Greedy: accept best unambiguous candidates first (a starting point
        for curation, not a replacement for it)."""
        done = []
        for c in self.candidates:
            tang = set(self.graph.tangles())
            if (c["score"] < min_score or "ambiguous" in c["flags"] or c["a"][0] == c["b"][0]
                    or c["a"][0] in tang or c["b"][0] in tang):
                continue
            try:
                self.project.join(c["a"], c["b"], c["gap"], reason=f"auto score={c['score']}",
                                  evidence=c["evidence"], via=c["via"])
                done.append(c["id"])
            except Exception:
                pass
        return {"joined": done}
