"""Load inputs, compute all evidence, answer editor queries."""
import os
from .gfa import AssemblyGraph
from .seqs import read_fasta
from .alignments import (AlignmentStore, read_alignments, read_paf, hits_from_alignments,
                         render_read, PLATFORM_DEFAULTS)
from . import softclip, longreads, evidence, circular
from .hic import HiC, read_pairs
from .project import Project


class Session:
    def __init__(self, gfa=None, fasta=None, alignments=(), paf=(), pairs=None, project=None,
                 hic_end_bin=20000, start_fasta=None):
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

        self.start_genes = read_fasta(start_fasta) if start_fasta else {}
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
        # snap read-based gap estimates to an exact sequence overlap where one exists
        for c in self.candidates:
            c["gap_reads"] = c["gap"]
            if c["score"] >= 2 and not c["via"]:
                g = circular.refine_gap(self.seqs, tuple(c["a"]), tuple(c["b"]), c["gap"])
                if g != c["gap"]:
                    c["gap"] = g
                    c["flags"].append(f"exact-overlap:{-g}")
        self.circular = circular.circularity_report(self.graph, self.candidates, self.seqs)
        # an exact duplicated end overlap beats the bucketed soft-clip gap estimate
        for r in self.circular:
            c = next((c for c in self.candidates if c["id"] == f"{r['contig']}:L~{r['contig']}:R"), None)
            if c:
                c["circular"] = True
                if "circular" not in c["flags"]:
                    c["flags"].append("circular")
                if r["end_overlap"]:
                    c["gap"] = -r["end_overlap"]
        self.cand_by_id = {c["id"]: c for c in self.candidates}
        self._explain_hotspots()

    def _explain_hotspots(self, tol=50):
        """Reads pile up (and clip) at the inner edge of a duplicated end
        overlap because they align to either copy; that's not a mis-join."""
        bounds = {}
        for r in self.circular:
            if r["end_overlap"]:
                L, o = r["length"], r["end_overlap"]
                bounds.setdefault(r["contig"], []).extend([o, L - o])
        for c in self.candidates:
            g = c["gap"]
            if g is not None and g < 0 and c["score"] >= 6:
                for name, side in (c["a"], c["b"]):
                    bounds.setdefault(name, []).append(-g if side == "L" else self.lengths[name] + g)
        for h in self.hotspots:
            if any(abs(h["pos"] - b) <= tol for b in bounds.get(h["contig"], [])):
                h["explained"] = "end-overlap boundary"
        self.hotspots.sort(key=lambda h: ("explained" in h, -h["reads"]))

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
            "circular": self.circular,
            "start_genes": list(self.start_genes),
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
            if "gap" in req and req["gap"] is not None:
                gap = int(req["gap"])
            elif c and c["via"]:
                gap = c["gap"]
            else:
                gap = circular.refine_gap(self.seqs, tuple(a), tuple(b), c["gap"] if c else 100)
            return {"scaffold": P.join(a, b, gap, reason=req.get("reason", ""),
                                       evidence=c["evidence"] if c else None,
                                       via=req.get("via", c["via"] if c else []))}
        if op == "break":
            return P.break_scaffold(req["scaffold"], int(req["after"]))
        if op == "split":
            return P.split_contig(req["piece"], int(req["pos"]), reason=req.get("reason", ""))
        if op == "circularise":
            return P.circularise(req["scaffold"], int(req.get("gap", 0)), reason=req.get("reason", ""))
        if op == "linearise":
            return P.linearise(req["scaffold"])
        if op == "set_origin":
            if "gene" in req:
                hit = self.locate_start(req["scaffold"], req["gene"])
                if not hit:
                    raise ValueError(f"{req['gene']} not found in {req['scaffold']}")
                P.set_origin(req["scaffold"], hit["pos"], hit["strand"], reason=f"{req['gene']} votes={hit['votes']}")
                return hit
            return P.set_origin(req["scaffold"], int(req["pos"]), req.get("strand", "+"), reason="manual")
        if op == "discard":
            return P.discard(req["scaffold"], reason=req.get("reason", ""))
        if op == "flip":
            return P.flip(req["scaffold"])
        if op == "undo":
            return P.undo()
        if op == "auto":
            return self.auto_scaffold(float(req.get("min_score", 6)), circularise=req.get("circularise", True))
        raise ValueError(f"unknown op {op}")

    def locate_start(self, scaffold, gene=None):
        s = self.project._scaffold(scaffold)
        seq = self.project.scaffold_seq(s, rotate=False)
        best = None
        for name, q in self.start_genes.items():
            if gene and name != gene:
                continue
            hit = circular.find_start(seq, q)
            if hit and (best is None or hit["votes"] > best["votes"]):
                best = {**hit, "gene": name}
        return best

    def auto_scaffold(self, min_score=6.0, circularise=True, min_circ_score=4.0):
        """Greedy: accept best unambiguous candidates first, then close circles
        and rotate to the start gene (a starting point for curation)."""
        done, rotated = [], []
        was_circular = {sc["name"] for sc in self.project.scaffolds if sc.get("circular")}
        tang = set(self.graph.tangles())
        for c in self.candidates:
            if c.get("circular"):
                continue
            if (c["score"] < min_score or "ambiguous" in c["flags"] or c["a"][0] == c["b"][0]
                    or c["a"][0] in tang or c["b"][0] in tang):
                continue
            try:
                gap = c["gap"] if c["via"] else circular.refine_gap(self.seqs, tuple(c["a"]), tuple(c["b"]), c["gap"])
                self.project.join(c["a"], c["b"], gap, reason=f"auto score={c['score']}",
                                  evidence=c["evidence"], via=c["via"])
                done.append(c["id"])
            except Exception:
                pass
        if circularise:
            # close any scaffold whose two terminal ends are joined by a candidate
            for c in self.candidates:
                if c["score"] < min_circ_score:
                    continue
                try:
                    sa, side_a = self.project._scaffold_end(tuple(c["a"]))
                    sb, side_b = self.project._scaffold_end(tuple(c["b"]))
                except Exception:
                    continue
                if sa != sb or side_a == side_b:
                    continue
                sc = self.project.scaffolds[sa]
                if sc.get("circular"):
                    continue
                # closing join: orient so the scaffold's R terminus comes first
                ea, eb = (c["a"], c["b"]) if side_a == "R" else (c["b"], c["a"])
                gap = circular.refine_gap(self.seqs, tuple(ea), tuple(eb), c["gap"])
                self.project.circularise(sc["name"], gap, reason=f"auto score={c['score']}")
        # scaffolds closed in either pass (a join between a scaffold's own two
        # termini circularises it)
        circles = [sc["name"] for sc in self.project.scaffolds
                   if sc.get("circular") and sc["name"] not in was_circular]
        if circularise:
            for name in circles:
                hit = self.locate_start(name) if self.start_genes else None
                if hit:
                    self.project.set_origin(name, hit["pos"], hit["strand"],
                                            reason=f"{hit['gene']} votes={hit['votes']}")
                    rotated.append({"scaffold": name, **hit})
        return {"joined": done, "circularised": circles, "rotated": rotated}
