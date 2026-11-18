# gapgraph — a Gap5-style finishing editor rebuilt around the assembly graph

Staden/Gap5 made finishing a *linear* job: one contig and its read pileup at a
time, with joins found by eye. Modern assemblies (hifiasm, Flye, Verkko) come out
as **graphs**, and the evidence for closing them comes from several sources at once:
long reads (ONT, PacBio CLR/HiFi), short reads, soft-clipped read ends and Hi-C.
gapgraph puts all of that evidence on one shared key, the **contig end**, and lets
you curate it in a single view.

| Gap5 idea | gapgraph equivalent |
|---|---|
| contig list | assembly graph view (Bandage-like), contig ends as nodes |
| "Find internal joins" | ranked **join candidates** that merge graph links, long-read bridges, soft-clip matches and Hi-C end contacts |
| contig editor | Gap5-style pileup with **soft-clipped bases shown past the contig end** (orange) and a **junction view** that shows both ends at once |
| break contig | **break hotspots**: piles of reads clipping at the same column inside a contig (chimeras, repeat boundaries) |
| repeats by hand | long-read **paths** through tangle nodes (X→R→Y). A join can carry a *copy* of the repeat |
| database + undo | JSON project file with an undo stack and audit log. Alignments stay in their own (indexed) files |
| output | FASTA + AGP 2.1 (gap evidence `proximity_ligation` for Hi-C-only joins) + edit log |

## Quick start (synthetic data with a known answer)

```bash
pip install -e .            # pure Python; `pip install -e .[bam]` adds pysam for BAM/CRAM
gapgraph simulate demo/
gapgraph serve --gfa demo/assembly.gfa --aln illumina:demo/illumina.sam \
  --aln hifi:demo/hifi.sam --aln ont:demo/ont.sam --pairs demo/hic.pairs --project demo/project.json
# open http://127.0.0.1:8765
```

The simulated genome is `A R B R C D E`. It is built to contain:

- a collapsed repeat `R`, which forms a tangle in the graph;
- a C|D gap with no graph link;
- contig E stored reverse-complemented;
- a 2 kb chimeric insert inside A.

The test suite checks that gapgraph finds all of these. After splitting the chimera and running auto-scaffold, the exported scaffold is **identical to the true genome**.

Keys: `j`/`k` next/prev join · `a` accept · `x` reject · `u` undo · `←`/`→` scroll · `[`/`]` jump to contig ends.

## Real data

```bash
# long reads -> contigs: keep soft-clips on supplementary alignments too (-Y)
minimap2 -ax map-ont -Y contigs.fa ont.fq.gz  | samtools sort -o ont.bam  && samtools index ont.bam
minimap2 -ax map-hifi -Y contigs.fa hifi.fq.gz | samtools sort -o hifi.bam && samtools index hifi.bam
bwa mem -Y contigs.fa R1.fq.gz R2.fq.gz | samtools sort -o ill.bam
# Hi-C -> .pairs (e.g. bwa mem -5SP | pairtools parse | pairtools dedup)

gapgraph report --gfa asm.gfa --fasta contigs.fa --aln ont:ont.bam --aln hifi:hifi.bam \
  --aln illumina:ill.bam --pairs hic.pairs
gapgraph serve  ...same inputs... --project asm.gapgraph.json --out asm.curated
```

Platforms (`illumina|hifi|clr|ont`) set the defaults for minimum clip length, k-mer size, end window and minimum MAPQ (`alignments.PLATFORM_DEFAULTS`). You can also pass long-read hits as `--paf ont:reads.paf`.

## Circular bacterial genomes (chromosome + plasmids)

A finished bacterial genome is a set of **closed circles**, each with no duplicated sequence, and each chromosome **starting at *dnaA* on the forward strand**. gapgraph handles this as follows:

- **Reads that wrap the origin count as evidence.** A long read that runs off a contig's R end and continues at its L end is recorded as a `contig:L ~ contig:R` candidate, flagged `circular`. The same applies to soft-clips at one end that match the other end. For a chromosome left in several contigs, the closing join is simply the candidate between the scaffold's two terminal ends.
- **Duplicated end overlaps are found exactly.** Canu, and sometimes Flye or hifiasm, leave the same sequence at both ends of a circular contig. Read-based gap estimates are snapped to the exact suffix/prefix overlap (flag `exact-overlap:N`), and closing the circle trims it.
- **Overlap edges are not reported as mis-joins.** Reads pile up and clip at the inner edge of a duplicated overlap. These hotspots are labelled `end-overlap boundary` instead of being shown as break candidates.
- **Start rotation.** `--start-fasta` takes one or more start genes, usually *dnaA* nucleotide sequences from the same or a related species; a few % divergence is fine. gapgraph finds the gene on either strand, including when it spans the current origin. It then rotates the circle, and reverse-complements it if needed, so the gene is base 1 and forward. The AGP is re-cut at the new origin.
- **Output.** FASTA headers carry `circular=true [topology=circular]`, which annotation tools such as Bakta read. The AGP marks circular objects with a comment, because AGP itself has no topology field.

```bash
# try it on a simulated chromosome (2 contigs, one reversed, overlapping
# ends, wrong start) + plasmid (duplicated 600 bp end overlap)
gapgraph simulate --circular demo_bact/
gapgraph serve --gfa demo_bact/assembly.gfa --aln illumina:demo_bact/illumina.sam \
  --aln ont:demo_bact/ont.sam --start-fasta demo_bact/dnaA_related.fasta
```

In the editor, open **Scaffolds**, which lists the circularity evidence per contig. Accept the joins, then use **Close circle**, then **start: at gene**. To do the same without the editor:

```bash
gapgraph auto --gfa demo_bact/assembly.gfa --aln illumina:demo_bact/illumina.sam \
  --aln ont:demo_bact/ont.sam --start-fasta demo_bact/dnaA_related.fasta --out demo_bact/finished
```

On the simulation, both replicons come out identical to the truth, base for base (`tests/test_gapgraph.py::test_bacterial_auto_finish_is_exact`).

**Real data**

The typical case is a Flye, hifiasm or Canu assembly from ONT or HiFi reads, plus Illumina reads for polishing.

```bash
minimap2 -ax map-ont -Y asm.fasta ont.fq.gz | samtools sort -o ont.bam && samtools index ont.bam
bwa index asm.fasta && bwa mem -Y asm.fasta R1.fq.gz R2.fq.gz | samtools sort -o ill.bam && samtools index ill.bam
gapgraph report --gfa assembly_graph.gfa --fasta asm.fasta --aln ont:ont.bam --aln illumina:ill.bam --start-fasta dnaA.fna
```

Things to watch with real bacterial data:

- Install pysam (`pip install -e ".[bam]"`) to read BAM. A bacterial genome fits comfortably in memory.
- Map reads to the same sequences gapgraph loads. Flye's `assembly_graph.gfa` uses `edge_N` segment names, while its `assembly.fasta` uses `contig_N` names. Either run gapgraph with `--fasta assembly.fasta` only (no `--gfa`), or extract the graph segments with `awk '/^S/{print ">"$2"\n"$3}' assembly_graph.gfa > edges.fa` and map the reads to `edges.fa`. Flye's `assembly_info.txt` also says which contigs it considers circular, which is a useful cross-check.
- Small plasmids can be shorter than the long reads, so a read may wrap the plasmid several times. These reads show up as repeated `plasmid+ -> plasmid+` paths, which is expected.
- Polish after rotation and closing, not before: polishers handle the new origin best when it is far from a contig end.

## How the evidence is scored (`evidence.py`)

- **graph link**: +3.
- **long-read bridges**: consecutive hits on one read leave contig X through one end and enter contig Y through another. Score +2·log2(1+n). Gap = read distance minus the unaligned contig ends. A read that leaves mid-contig is counted as a conflict, not as support.
- **soft-clips**: each clip is oriented outward and k-mer matched against the inward sequence of every other end. A consistent diagonal gives the partner and the gap. Score +log2(1+n).
- **Hi-C**: contacts are counted per contig-end pair. Score +6·share·(0.5+orientation fraction), where:
  - share = the geometric mean of the fraction of each end's contacts that go to the other end;
  - orientation fraction = this end pair's share of the 4 possible end pairs for these two contigs.
- **flags**:
  - `ambiguous`: another candidate on the same end scores within 60%;
  - `via:R+`: the join goes through a repeat;
  - `outcompeted`: there is a graph link, but the data prefers another partner;
  - `unsupported-link`;
  - `hic-only`;
  - `gap-to-fill`.

## Limits of this prototype

- Alignments are held in memory per contig. For large genomes, the next step is pysam region fetches from indexed BAM.
- The soft-clip consensus is a naive column vote. It works for HiFi and Illumina, but is only indicative for ONT/CLR; a POA consensus would be the next step.
- Gap filling only writes Ns. Splicing in the clip or read consensus is not yet implemented.
- There is no base-level editing. Polishing belongs to dedicated tools.

## References

- Bonfield JK, Whitwham A (2010) Gap5—editing the billion fragment sequence assembly. *Bioinformatics* 26:1699–1703. doi:10.1093/bioinformatics/btq268
- Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094–3100. doi:10.1093/bioinformatics/bty191
- Wick RR et al. (2015) Bandage: interactive visualization of de novo genome assemblies. *Bioinformatics* 31:3350–3352. doi:10.1093/bioinformatics/btv383
- Cheng H et al. (2021) Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nat Methods* 18:170–175. doi:10.1038/s41592-020-01056-5
- Hunt M et al. (2015) Circlator: automated circularization of genome assemblies using long sequencing reads. *Genome Biol* 16:294. doi:10.1186/s13059-015-0849-0 (an established standalone tool for the same closure/rotation step)
- Kolmogorov M et al. (2019) Assembly of long, error-prone reads using repeat graphs (Flye). *Nat Biotechnol* 37:540–546. PMID 30936562
