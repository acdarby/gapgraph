"""Load inputs, compute all evidence, answer editor queries."""
import os
from .gfa import AssemblyGraph
from .seqs import read_fasta, revcomp
from .alignments import (AlignmentStore, read_alignments, read_paf, hits_from_alignments,
                         render_read, PLATFORM_DEFAULTS)
from . import softclip, longreads, evidence, circular, pairs, joinedit
from .consensus import Pileup, DEFAULT_THRESHOLDS, PROBLEM_LABELS
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
        # read pairs (Gap5 "Find Read Pairs")
        self.insert_stats = {}
        self.pair_joins = {}
        for plat, st in self.stores.items():
            stats = pairs.insert_stats(st)
            if stats:
                self.insert_stats.update(stats)
                for k, v in pairs.spanning_joins(st, self.lengths, stats).items():
                    m = self.pair_joins.setdefault(k, {"pairs": 0, "gap": v["gap"]})
                    m["pairs"] += v["pairs"]
        self.other_contigs = pairs.other_contig_index(self.stores)
        self._pileups = {}
        self.candidates = evidence.build_candidates(self.graph, clip_joins, self.lr_joins, self.hic,
                                                    self.lr_paths, self.pair_joins)
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
            "insert_stats": self.insert_stats,
            "problem_types": PROBLEM_LABELS,
            "problem_defaults": DEFAULT_THRESHOLDS,
        }

    def region(self, contig, start, end, max_rows=400, overhang=20000):
        """Editor window.  The window may extend past either contig end so the
        soft-clipped bases of reads running off the end stay visible."""
        L = self.lengths[contig]
        start, end = max(-overhang, start), min(L + overhang, end)
        if end <= start:
            end = start + 1
        rows, row_end, reads = [], [], []
        for plat in sorted(self.stores):
            reads.extend(self.stores[plat].fetch(contig, start, end))
        reads.sort(key=lambda a: a.pos - a.lclip_soft)
        for a in reads:
            r = render_read(a, start, end)
            if not r["s"]:
                continue
            r["st"] = pairs.classify(a, self.insert_stats, self.other_contigs)
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
        cons, qual = self.consensus_window(contig, start, end)
        return {"contig": contig, "start": start, "end": end, "length": L,
                "ref": ref, "rows": rows, "depth": depth, "consensus": cons, "qual": qual}

    # ------------------------------------------------- consensus & problems
    def pileup(self, contig):
        P = self._pileups.get(contig)
        if P is None:
            alns = (a for st in self.stores.values() for a in st.by_ref.get(contig, []))
            P = self._pileups[contig] = Pileup(contig, self.lengths[contig]).add_alignments(alns)
        return P

    def consensus_window(self, contig, start, end):
        P = self.pileup(contig)
        L = self.lengths[contig]
        lo, hi = max(0, start), min(L, end)
        pad_l, pad_r = max(0, -start) if start < 0 else 0, max(0, end - L)
        cons = " " * pad_l + P.consensus_string(lo, hi) + " " * pad_r
        qual = [0] * pad_l + P.qual[lo:hi].tolist() + [0] * pad_r
        return cons, qual

    def problems(self, contig, kind, thr=None):
        thr = DEFAULT_THRESHOLDS[kind] if thr is None else float(thr)
        P = self.pileup(contig)
        return P.regions(P.masks(self.seqs[contig], kind, thr))

    def next_problem(self, contig, pos, kind, thr=None, direction=1):
        regs = self.problems(contig, kind, thr)
        if direction > 0:
            for s, e in regs:
                if s > pos:
                    return {"start": s, "end": e}
        else:
            for s, e in reversed(regs):
                if s < pos:
                    return {"start": s, "end": e}
        return None

    def problems_summary(self, kind, thr=None, limit=500):
        out = []
        for c in self.seqs:
            for s, e in self.problems(c, kind, thr):
                out.append({"contig": c, "start": s, "end": e})
                if len(out) >= limit:
                    return out
        return out

    # ------------------------------------------------- template display
    def templates(self, contig, start, end):
        res = pairs.templates(self.stores, contig, start, end, self.insert_stats, self.other_contigs)
        res["depth"] = pairs.depth_track(self.stores, contig, self.lengths[contig], start, end)
        res["length"] = self.lengths[contig]
        res["insert_stats"] = self.insert_stats
        return res

    # ------------------------------------------------- join editor
    def _join_ends(self, req):
        c = self.cand_by_id.get(req.get("id"))
        if c:
            return tuple(c["a"]), tuple(c["b"]), int(req.get("gap", c["gap"])), c
        return tuple(req["a"]), tuple(req["b"]), int(req.get("gap", 0)), None

    def joinview(self, req, start, end, max_rows=200):
        a, b, gap, c = self._join_ends(req)
        la, lb = self.lengths[a[0]], self.lengths[b[0]]
        a_rev, b_rev, a_off, b_off = joinedit.layout(a, b, gap, la, lb)
        out = {"a": list(a), "b": list(b), "gap": gap, "start": start, "end": end,
               "a_range": [a_off, 0], "b_range": [b_off, b_off + lb]}
        for key, (ctg, L, rev, off) in (("A", (a[0], la, a_rev, a_off)), ("B", (b[0], lb, b_rev, b_off))):
            # join window -> oriented contig cols -> original cols
            o0, o1 = start - off, end - off
            if rev:
                f0, f1 = L - o1, L - o0
            else:
                f0, f1 = o0, o1
            reg = self.region(ctg, f0, f1, max_rows=max_rows)
            rows = []
            for row in reg["rows"]:
                nr = []
                for r in row:
                    r2 = joinedit.orient_read(r, L, rev)
                    r2 = dict(r2, x=r2["x"] + off, ins=[[p + off, s] for p, s in r2["ins"]])
                    nr.append(r2)
                rows.append(nr)
            ref, cons, qual, depth = reg["ref"], reg["consensus"], reg["qual"], reg["depth"]
            if rev:
                ref = ref.translate(joinedit._COMP)[::-1]
                cons = cons.translate(joinedit._COMP)[::-1]
                qual, depth = qual[::-1], depth[::-1]
            out[key] = {"contig": ctg, "rev": rev, "rows": rows, "ref": ref, "consensus": cons,
                        "qual": qual, "depth": depth}
        return out

    def oriented_seq(self, end, right):
        s = self.seqs[end[0]]
        if right:   # end on the right
            return s if end[1] == "R" else revcomp(s)
        return s if end[1] == "L" else revcomp(s)

    def joinalign(self, req):
        a, b, gap, c = self._join_ends(req)
        res = joinedit.align_overlap(self.oriented_seq(a, True), self.oriented_seq(b, False), gap)
        res.update(a=list(a), b=list(b), gap=gap)
        return res

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
