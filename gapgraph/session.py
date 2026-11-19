"""Load inputs, compute all evidence, answer editor queries."""
import os
from .gfa import AssemblyGraph
from .seqs import read_fasta, revcomp
from .alignments import (AlignmentStore, read_alignments, read_paf, hits_from_alignments,
                         render_read, PLATFORM_DEFAULTS)
from . import softclip, longreads, evidence, circular, pairs, joinedit
from .consensus import Pileup, DEFAULT_THRESHOLDS, PROBLEM_LABELS
from .hic import HiC, read_pairs
from .project import Project, ProjectError
from .working import Working
from . import realign
from collections import Counter


class Session:
    def __init__(self, gfa=None, fasta=None, alignments=(), paf=(), pairs=None, project=None,
                 hic_end_bin=20000, start_fasta=None):
        """alignments / paf: lists of (path, platform) with platform in
        illumina | hifi | clr | ont."""
        if gfa:
            graph = AssemblyGraph.from_gfa(gfa)
            if fasta:
                for n, sq in read_fasta(fasta).items():
                    if n in graph.segments:
                        graph.segments[n].seq = sq
                        graph.segments[n].length = len(sq)
        else:
            graph = AssemblyGraph.from_fasta(read_fasta(fasta))
        seqs = {n: sg.seq for n, sg in graph.segments.items()}
        missing = [n for n, sq in seqs.items() if not sq]
        if missing:
            raise ValueError(f"no sequence for segments {missing[:5]} -- pass --fasta")
        # the inputs, never modified
        self.orig_graph, self.orig_seqs = graph, seqs
        self.orig_stores = {}
        for path, plat in alignments:
            self.orig_stores.setdefault(plat, AlignmentStore()).add(read_alignments(path, plat))
        self.orig_hits = []
        for path, plat in paf:
            self.orig_hits.extend(read_paf(path, plat))
        self.hic_raw = list(read_pairs(pairs)) if pairs else []
        self.hic_end_bin = hic_end_bin

        self.start_genes = read_fasta(start_fasta) if start_fasta else {}
        self.project = Project(dict(seqs), path=project)
        if project and os.path.exists(project):
            self.project.load()
        self._wcache = {}
        self.rebuild()

    def rebuild(self):
        """Derive the working assembly from the project and recompute evidence."""
        if len(self._wcache) > 2000:
            self._wcache.clear()
        W = Working(self.project, self.orig_graph, self.orig_seqs, self.orig_stores,
                    self.orig_hits, self.hic_raw, cache=self._wcache)
        self.W = W
        self.graph, self.seqs, self.lengths, self.stores = W.graph, W.seqs, W.lengths, W.stores
        self.hits = list(W.hits)
        for plat, st in self.stores.items():
            if plat != "illumina":
                self.hits.extend(hits_from_alignments(st.all(), self.lengths))
        self.hic = HiC(self.lengths, end_bin=self.hic_end_bin).add(W.hic_pairs) if self.hic_raw else None
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
        self._ext = {}
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
        for c in self.candidates:
            ca, cb = c["a"][0], c["b"][0]
            if ca in self.W.circular or cb in self.W.circular:
                if ca == cb:
                    c["done"] = True       # this circle is already closed
                elif "circular-contig" not in c["flags"]:
                    c["flags"].append("circular-contig")
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
            "working": {w: {"scaffold": self.W.scaffold_of[w], "circular": w in self.W.circular,
                            "parts": len(self.W.comps[w]),
                            "junctions": [off for off, c in self.W.comps[w][1:]]} for w in self.seqs},
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

    # ------------------------------------------------- clips beyond contig ends
    def _end_reads(self, w, side, tol=30):
        L, seq = self.lengths[w], self.seqs[w]
        lo, hi = (max(0, L - tol - 1), L) if side == "R" else (0, tol + 1)
        for st in self.stores.values():
            for a in st.fetch(w, lo, hi):
                if not a.seq:
                    continue
                if side == "R" and a.rclip_soft >= 3 and a.end >= L - tol:
                    yield a, a.seq[len(a.seq) - a.rclip_soft:], seq[a.end:L]
                elif side == "L" and a.lclip_soft >= 3 and a.pos <= tol:
                    yield a, revcomp(a.seq[:a.lclip_soft]), revcomp(seq[:a.pos])

    @staticmethod
    def _rkey(a):
        return (a.qname, a.pos, a.reverse, a.supplementary)

    def end_ext(self, w, side, beyond=None):
        """Align every soft clip at one contig end to the sequence beyond it.

        `beyond` (outward orientation) is what is believed to lie past the end:
        the partner contig in the Join Editor, the contig's own start for a
        circle, and otherwise the consensus of the clips themselves (Gap5's
        'hidden data' extension)."""
        key = (w, side, hash(beyond) if beyond is not None else None)
        if key in self._ext:
            return self._ext[key]
        reads = list(self._end_reads(w, side))
        seq = self.seqs[w]
        if beyond is None:
            if w in self.W.circular:
                beyond = seq[:60000] if side == "R" else revcomp(seq[-60000:])
            else:
                rems = []
                for a, clip, tail in reads:
                    q = realign.prefix_align(clip, tail)["q"] if tail else 0
                    rems.append(clip[q:])
                beyond = realign.clip_consensus(rems) if len(rems) >= 3 else max(rems, key=len, default="")
        out, votes = {}, []
        for a, clip, tail in reads:
            tgt = tail + beyond[:len(clip) + 200]
            al = realign.prefix_align(clip, tgt)
            chars, ins = realign.render_cols(clip[:al["q"]], al["cols"])
            out[self._rkey(a)] = {"chars": chars, "ins": ins, "raw": clip[al["q"]:], "ident": al["identity"]}
            for j, ch in enumerate(chars):
                k = j - len(tail)
                if k >= 0:
                    while len(votes) <= k:
                        votes.append(Counter())
                    votes[k][ch.upper()] += 1
        cons = []
        for v in votes:
            n = sum(v.values())
            cons.append(v.most_common(1)[0][0].lower() if n >= 2 else " ")
        res = {"reads": out, "cons": "".join(cons).rstrip(), "beyond": beyond}
        self._ext[key] = res
        return res

    def _render(self, a, start, end, ext_r, ext_l):
        """render_read with end clips replaced by their aligned versions."""
        base = render_read(a, -10 ** 9, 10 ** 9)
        x0, s, ins, unal = base["x"], base["s"], list(base["ins"]), []
        _comp = str.maketrans("ACGTNacgtn*", "TGCANtgcan*")
        if ext_r:
            s = s[:len(s) - a.rclip_soft] if a.rclip_soft else s
            ch = ext_r["chars"].lower()
            ins += [[a.end + j, b.lower()] for j, b in ext_r["ins"]]
            if ext_r["raw"]:
                unal.append([x0 + len(s) + len(ch), x0 + len(s) + len(ch) + len(ext_r["raw"])])
            s = s + ch + ext_r["raw"].lower()
        if ext_l:
            if a.lclip_soft:
                s = s[a.lclip_soft:]
            ch = ext_l["chars"].translate(_comp)[::-1].lower()
            raw = revcomp(ext_l["raw"]).lower()
            ins += [[a.pos - j, revcomp(b).lower()] for j, b in ext_l["ins"]]
            s = raw + ch + s
            x0 = a.pos - len(ch) - len(raw)
            if raw:
                unal.append([x0, x0 + len(raw)])
        lo, hi = max(start, x0), min(end, x0 + len(s))
        r = dict(base)
        r.update(x=lo, s=s[lo - x0:hi - x0] if hi > lo else "", ins=[i for i in ins if lo <= i[0] < hi],
                 u=[[max(lo, u0), min(hi, u1)] for u0, u1 in unal if u1 > lo and u0 < hi])
        return r

    def region(self, contig, start, end, max_rows=400, overhang=20000, beyond=None):
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
        beyond = beyond or {}
        ext = {}
        if end > L - 60:
            ext["R"] = self.end_ext(contig, "R", beyond.get("R"))
        if start < 60:
            ext["L"] = self.end_ext(contig, "L", beyond.get("L"))
        for a in reads:
            k = self._rkey(a)
            er = ext["R"]["reads"].get(k) if "R" in ext else None
            el = ext["L"]["reads"].get(k) if "L" in ext else None
            r = self._render(a, start, end, er, el) if (er or el) else render_read(a, start, end)
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
        # extend the consensus line past the contig ends with the clip consensus
        cl = list(cons)
        if "R" in ext:
            for k, ch in enumerate(ext["R"]["cons"]):
                i = L + k - start
                if 0 <= i < len(cl):
                    cl[i] = ch
                    qual[i] = 75
        if "L" in ext:
            for k, ch in enumerate(ext["L"]["cons"]):
                i = -1 - k - start
                if 0 <= i < len(cl) and ch != " ":
                    cl[i] = ch.translate(str.maketrans("acgtn*", "tgcan*"))
                    qual[i] = 75
        cons = "".join(cl)
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
        A_or, B_or = self.oriented_seq(a, True), self.oriented_seq(b, False)
        cap = 60000
        if gap <= 0:
            beyond_a = B_or[-gap:][:cap]
            beyond_b = revcomp(A_or[:la + gap][-cap:])
        else:
            fill_a = self.end_ext(a[0], a[1])["beyond"][:gap].upper()
            fill_b = self.end_ext(b[0], b[1])["beyond"][:gap].upper()
            beyond_a = (fill_a + "N" * (gap - len(fill_a)) + B_or)[:cap]
            beyond_b = (fill_b + "N" * (gap - len(fill_b)) + revcomp(A_or))[:cap]
        bey = {"A": {a[1]: beyond_a}, "B": {b[1]: beyond_b}}
        for key, (ctg, L, rev, off) in (("A", (a[0], la, a_rev, a_off)), ("B", (b[0], lb, b_rev, b_off))):
            # join window -> oriented contig cols -> original cols
            o0, o1 = start - off, end - off
            if rev:
                f0, f1 = L - o1, L - o0
            else:
                f0, f1 = o0, o1
            reg = self.region(ctg, f0, f1, max_rows=max_rows, beyond=bey[key])
            rows = []
            for row in reg["rows"]:
                nr = []
                for r in row:
                    r2 = joinedit.orient_read(r, L, rev)
                    u = r2.get("u", [])
                    if rev:
                        u = [[L - u1, L - u0] for u0, u1 in u]
                    r2 = dict(r2, x=r2["x"] + off, ins=[[p + off, s] for p, s in r2["ins"]],
                              u=[[u0 + off, u1 + off] for u0, u1 in u])
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
        return self.hic.matrix(list(self.seqs))

    def _project_via(self, via):
        out = []
        for v in via or ():
            w, o = v[:-1], v[-1]
            s = self.project._scaffold(self.W.scaffold_of[w])
            if len(s["parts"]) != 1:
                return []
            part = s["parts"][0]
            pc = self.project.pieces[part["piece"]]
            if pc["start"] != 0 or pc["end"] != len(self.project.seqs[pc["src"]]):
                return []
            out.append(pc["src"] + (o if part["orient"] == "+" else ("-" if o == "+" else "+")))
        return out

    def op(self, req):
        before = set(self.seqs)
        res = self._op(req)
        self.rebuild()
        return {"result": res, "new": sorted(set(self.seqs) - before), "gone": sorted(before - set(self.seqs))}

    def _op(self, req):
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
            via = self._project_via(req.get("via", c["via"] if c else []))
            pa, pb = self.W.to_project_end(tuple(a)), self.W.to_project_end(tuple(b))
            return {"scaffold": P.join(pa, pb, gap, reason=req.get("reason", ""),
                                       evidence=c["evidence"] if c else None, via=via)}
        if op == "break":
            return P.break_scaffold(req["scaffold"], int(req["after"]))
        if op == "split":
            if "piece" in req:
                return P.split_contig(req["piece"], int(req["pos"]), reason=req.get("reason", ""))
            w, x = req["contig"], int(req["pos"])
            if w in self.W.circular:
                raise ProjectError(f"{w} is circular: linearise it before splitting")
            loc = self.W.locate(w, x)
            if loc is None:
                raise ProjectError("split position outside the contig")
            i, piece, rel, off = loc
            if rel == 0 and i > 0:     # exactly on a join: break the scaffold there
                return P.break_scaffold(self.W.scaffold_of[w], i - 1)
            return P.split_contig(piece, rel, reason=req.get("reason", ""))
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
            return self._auto(float(req.get("min_score", 6)), circularise=req.get("circularise", True))
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
        res = self._auto(min_score, circularise, min_circ_score)
        self.rebuild()
        return res

    def _auto(self, min_score=6.0, circularise=True, min_circ_score=4.0):
        """Greedy: accept best unambiguous candidates first, then close circles
        and rotate to the start gene (a starting point for curation)."""
        done, rotated = [], []
        # translate every candidate to project ends *before* the layout changes
        tr = {}
        for c in self.candidates:
            if c.get("done"):
                continue
            try:
                tr[c["id"]] = (self.W.to_project_end(tuple(c["a"])), self.W.to_project_end(tuple(c["b"])),
                               self._project_via(c["via"]))
            except Exception:
                pass
        was_circular = {sc["name"] for sc in self.project.scaffolds if sc.get("circular")}
        tang = set(self.graph.tangles())
        for c in self.candidates:
            if c.get("circular") or c["id"] not in tr:
                continue
            if (c["score"] < min_score or "ambiguous" in c["flags"] or c["a"][0] == c["b"][0]
                    or c["a"][0] in tang or c["b"][0] in tang):
                continue
            try:
                gap = c["gap"] if c["via"] else circular.refine_gap(self.seqs, tuple(c["a"]), tuple(c["b"]), c["gap"])
                pa, pb, via = tr[c["id"]]
                self.project.join(pa, pb, gap, reason=f"auto score={c['score']}",
                                  evidence=c["evidence"], via=via)
                done.append(c["id"])
            except Exception:
                pass
        if circularise:
            # close any scaffold whose two terminal ends are joined by a candidate
            for c in self.candidates:
                if c["score"] < min_circ_score or c["id"] not in tr:
                    continue
                pa, pb, _ = tr[c["id"]]
                try:
                    sa, side_a = self.project._scaffold_end(pa)
                    sb, side_b = self.project._scaffold_end(pb)
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
