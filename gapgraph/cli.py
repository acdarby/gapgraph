"""Command line: gapgraph {simulate,report,serve}."""
import argparse
import json
import sys


def _session(a):
    from .session import Session
    aln = [tuple(x.split(":", 1)[::-1]) if ":" in x else (x, "ont") for x in (a.aln or [])]
    paf = [tuple(x.split(":", 1)[::-1]) if ":" in x else (x, "ont") for x in (a.paf or [])]
    return Session(gfa=a.gfa, fasta=a.fasta, alignments=aln, paf=paf, pairs=a.pairs,
                   project=a.project, hic_end_bin=a.hic_end_bin, start_fasta=a.start_fasta,
                   foldback_mode=a.foldback_mode)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="gapgraph", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("simulate", help="write a synthetic finishing problem with known truth")
    s.add_argument("outdir")
    s.add_argument("--scale", type=float, default=0.5)
    s.add_argument("--seed", type=int, default=7)
    s.add_argument("--circular", action="store_true", help="bacterial chromosome + plasmid scenario")
    s.add_argument("--foldbacks", action="store_true",
                   help="(with --circular) add 5%% ONT foldback reads and a hairpin artefact on chr_1's R end")

    for name, hlp in (("report", "print ranked join candidates + break hotspots"),
                      ("serve", "start the interactive editor"),
                      ("auto", "greedy auto-scaffold then export")):
        p = sub.add_parser(name, help=hlp)
        p.add_argument("--gfa", help="assembly graph (GFA1; hifiasm/Flye/Verkko/SPAdes)")
        p.add_argument("--fasta", help="contig FASTA (needed if GFA has no sequences)")
        p.add_argument("--aln", action="append",
                       help="PLATFORM:file.sam|bam  (PLATFORM = illumina|hifi|clr|ont); repeatable")
        p.add_argument("--paf", action="append", help="PLATFORM:file.paf long-read hits; repeatable")
        p.add_argument("--pairs", help="Hi-C .pairs on contig coordinates")
        p.add_argument("--project", default="gapgraph_project.json", help="edit state (created if absent)")
        p.add_argument("--hic-end-bin", type=int, default=20000)
        p.add_argument("--foldback-mode", choices=("trim", "keep", "drop"), default="trim",
                       help="foldback reads: trim the folded arm (default), keep, or drop the read")
        p.add_argument("--start-fasta", help="start gene(s), e.g. dnaA nucleotide seqs from related strains, "
                                             "used to rotate circular replicons")
        p.add_argument("--out", default="gapgraph_out", help="export prefix")
        if name == "serve":
            p.add_argument("--port", type=int, default=8765)
            p.add_argument("--host", default="127.0.0.1")
        if name == "auto":
            p.add_argument("--min-score", type=float, default=6.0)
            p.add_argument("--no-circularise", action="store_true")
    a = ap.parse_args(argv)

    if a.cmd == "simulate":
        d = a.outdir
        if a.circular:
            from .simulate import simulate_circular
            r = simulate_circular(a.outdir, seed=a.seed, scale=a.scale,
                                  foldback_reads=0.05 if a.foldbacks else 0.0, foldback_contig=a.foldbacks)
            print(json.dumps(r["truth"], indent=1))
            print(f"\ntry:\n  gapgraph serve --gfa {d}/assembly.gfa --aln illumina:{d}/illumina.sam "
                  f"--aln ont:{d}/ont.sam --start-fasta {d}/dnaA_related.fasta --project {d}/project.json")
            return
        from .simulate import simulate
        r = simulate(a.outdir, seed=a.seed, scale=a.scale)
        print(json.dumps(r["truth"], indent=1))
        print(f"\ntry:\n  gapgraph serve --gfa {d}/assembly.gfa --aln illumina:{d}/illumina.sam "
              f"--aln hifi:{d}/hifi.sam --aln ont:{d}/ont.sam --pairs {d}/hic.pairs --project {d}/project.json")
        return
    if not (a.gfa or a.fasta):
        sys.exit("need --gfa and/or --fasta")
    S = _session(a)
    from . import foldback
    if a.cmd == "report":
        print("rank\tscore\tjoin\tevidence\tgap\tlr\tclip\thic\tflags")
        for i, c in enumerate(S.candidates[:50], 1):
            print(f"{i}\t{c['score']}\t{c['id']}\t{','.join(c['evidence'])}\t{c['gap']}\t"
                  f"{c['longread_reads']}\t{c['softclip_reads']}\t{c['hic_contacts']}\t{','.join(c['flags'])}")
        print("\n# break hotspots (internal soft-clip piles)")
        for h in S.hotspots[:30]:
            print(f"{h['contig']}\t{h['pos']}\t{h['clip_side']}\t{h['reads']} reads\t{h['platform']}\t"
                  f"{h.get('explained', '')}")
        print("\n# foldback reads (platform, reads, foldbacks, fraction, split-alignment, folded-clip, median arm,"
              " not counted: inside a contig foldback arm / reference palindrome piles)")
        for plat, ev in S.fb_events.items():
            sm = foldback.summarise(ev, S.raw_stores[plat], S.fb_stats[plat])
            print(f"{plat}\t{sm['reads']}\t{sm['foldbacks']}\t{sm['fraction']}\t{sm['split']}\t{sm['clip']}\t"
                  f"{sm['median_arm']}\t{sm['at_contig_foldback']}/{sm['reference_palindrome_piles']}")
        print("\n# foldback contig ends (contig, side, fold point, arm, identity, reads spanning fold, depth, recommendation)")
        for f in S.contig_foldbacks:
            print(f"{f['contig']}\t{f['side']}\t{f['fold']}\t{f['arm']}\t{f['identity']}\t"
                  f"{f['spanning_reads']}\t{f['depth']}\t{f['recommend']}")
        print("\n# circularity (contig, length, graph self-link, end overlap, wrapping long reads, wrapping clips, score)")
        for r in S.circular:
            print(f"{r['contig']}\t{r['length']}\t{r['graph_self_link']}\t{r['end_overlap']}\t"
                  f"{r['wrap_longreads']}\t{r['wrap_softclips']}\t{r['score']}\t{r['name_hint'] or ''}")
        print("\n# long-read paths through >=3 contigs")
        for p in S.lr_paths[:20]:
            print(f"{p['reads']}\t{' -> '.join(p['path'])}")
    elif a.cmd == "auto":
        print(json.dumps(S.auto_scaffold(a.min_score, circularise=not a.no_circularise), indent=1))
        print("\n".join(S.project.export(a.out)))
    else:
        from .server import serve
        serve(S, a.host, a.port, a.out)


if __name__ == "__main__":
    main()
