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
| Join Editor | one locked canvas: contig A's reads above its consensus, a **differences line**, contig B's consensus and reads below; Align + overlap length / % identity before committing |
| consensus + "Search by ..." | quality-weighted consensus with per-column confidence; **Ctrl-S / Ctrl-R** jump to the next problem: low consensus quality, contig ≠ reads, excess discrepancies, second allele, low / high coverage |
| Template Display / Find Read Pairs | whole-contig plot of every template at height = insert (or read) length, coloured by status (grey by MAPQ, **blue** unpaired, **orange** mate in another contig, **red** bad orientation, **green** wraps a circular end) + read/template depth; spanning pairs become join evidence |
| output | FASTA + AGP 2.1 (gap evidence `proximity_ligation` for Hi-C-only joins) + edit log |

## Quick start (synthetic data with a known answer)

```bash
pip install -e .            # needs numpy + edlib; `pip install -e ".[bam]"` adds pysam for BAM/CRAM
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

## Edits change the data (as in Gap5)

**Splits and joins create real contigs.** gapgraph keeps your BAM/SAM files read-only. After every operation (join, break, flip, close circle, set start, undo) it rebuilds a *working assembly* from the project layout:

- **One contig per scaffold.** Each scaffold becomes a single contig, with overlaps trimmed and circles rotated. The graph, editor, template display, problem search and join candidates all switch to the new contigs straight away.
- **Reads follow the edit.** Every input alignment is cut to the piece it falls in. Bases outside the piece become soft clips, just as an aligner would report them against the new contig. Alignments are reverse-complemented when a piece is used on the minus strand, and mates are lifted too, so pairs that spanned two contigs become consistent pairs after the join.
- **Break Contig** (✂ in the editor, or *Break* on a hotspot): reads that crossed the break become soft-clipped on both sides. The evidence for re-joining is then visible straight away; if the break was wrong, the top candidate puts it back.
- **Join:** soft clips that ran off the old contig ends are re-aligned (edlib) into the sequence on the other side of the junction. Reads therefore span the join and the consensus is continuous across it. A read aligned to both old contigs (primary and supplementary) is kept once.

**Soft clips are aligned, not just shown.** Clipped bases are Gap5's "hidden data". They are re-aligned to whatever lies beyond the contig end:

- **Contig editor:** they align to the consensus of all the clips at that end (iterated star alignment). That consensus extends the consensus line past the contig end in lower case.
- **Circular contig:** they align to the contig's own start.
- **Join Editor:** they align to the partner contig at the current gap. Contig A's clips sit in register under contig B's sequence and vice versa, and mismatches are highlighted. Changing the gap re-aligns them.
- **Unaligned tails:** clip sequence past the point where the alignment stops being supported (X-drop, typically adapter or junk) is drawn faded.

## Foldbacks (reads and contig ends)

The editor's **Breaks** tab, `gapgraph report` and auto-finish all handle foldbacks.

**Foldback reads.** These are mostly a nanopore artefact: the complementary strand of the same molecule follows the template through the pore, so the read is a sequence, a short hairpin, then the reverse complement. Aligners report it in one of two ways:
- a primary and supplementary pair on **opposite strands of the same place**;
- one alignment whose soft clip is the reverse complement of the reference just before the clip point.

Left in, the folded arm double-counts depth and looks like evidence for inversions, mis-joins and false joins. gapgraph detects both forms on load and handles them with `--foldback-mode` (also switchable in the editor):

| Mode | Effect |
|---|---|
| `trim` (default) | keep the longer arm, drop the folded arm's alignment and hard-clip its bases |
| `keep` | leave everything as it is, for inspection |
| `drop` | remove the read entirely |

Two cases are deliberately *not* counted as foldback reads:
- reads inside the arm of a contig foldback, which only look folded because the reference is;
- piles of folded reads at one position, which point to a palindrome in the reference. Real foldback reads fall at random positions.

**Foldback contig ends.** Assemblers sometimes finish a contig with a hairpin, where the last n bases are the reverse complement of the n before them.
- **Detection:** a fold at point f makes the window match its own reverse complement along a single diagonal (S[i] = R[i + W − 2f]). One k-mer pass therefore finds the fold point and arm length, including when there is a short loop between the arms, and an alignment then confirms the inversion.
- **Read evidence:** each call counts reads that run straight through the fold.
  - An assembler hairpin has none, so it is marked **remove**.
  - Reads spanning the fold mean it may be real, e.g. a hairpin telomere on a linear replicon such as a *Borrelia* linear plasmid, or a terminal inverted repeat. It is then marked **check** and left alone by auto-finish.
- **Removal:** **Remove arm** cuts the contig at the fold and discards the arm in a single edit that one undo reverses. Auto-finish removes confident hairpins before joining.

```bash
gapgraph simulate --circular --foldbacks demo_fb/     # 5% ONT foldback reads + a 1.5 kb hairpin on chr_1
gapgraph report --gfa demo_fb/assembly.gfa --aln illumina:demo_fb/illumina.sam --aln ont:demo_fb/ont.sam
```

On that simulation, gapgraph finds exactly the simulated foldback reads and the hairpin's fold point and arm length. Auto-finish still reproduces both replicons base for base.

## The Gap5 tools, rebuilt

These follow the Gap5 manual (Staden package Texinfo sources, `staden/trunk/doc/manual/gap5*.texi` in the Staden SourceForge repository), adapted to graph assemblies and mixed short- and long-read data.

**Join Editor** (select a join; `a` or **Join…**)
- In Gap5 the Join Editor is "a pair of contig editor displays stacked above one another… flipped in Y so that the consensus appears at the bottom… separated only by a differences line".
- In gapgraph, contig A is oriented so its joining end sits at x = 0 and contig B starts at x = gap. Reads are reverse-complemented server-side when a contig is used on its minus strand, so both halves always scroll in register.
- **Align** runs a banded alignment of the overlapping consensus sequences and prints a Gap5-style alignment.
- **Join…** reports the overlap length, % identity, mismatches and indels before you commit.
- Editing the gap box re-registers contig B live.

**Consensus and problem search** (editor bar, or the **Problems** tab for all contigs)
- Each read base adds evidence `ln((1−e)·3/e)` for its base (or `*` for a deletion), where e comes from the base quality, or from a per-platform default when there is none.
- The consensus is the best-supported symbol. Its confidence is Q = −10·log10(1 − P(best)), capped at 100. On the consensus line, a darker background means lower confidence.
- Search types:

  | Type | Finds |
  |---|---|
  | `consq` | Q below the threshold (Gap5's default search, threshold 30) |
  | `contig` | the contig base disagrees with the read consensus, including indels: unpolished errors |
  | `discrep` | more disagreeing reads than the base qualities predict |
  | `het` | a second allele, e.g. a mixed population or a collapsed repeat copy |
  | `lowcov` | depth below the threshold |
  | `highcov` | depth above N × the median. Collapsed rRNA operons and IS elements show up here; the manual notes "regions of extreme depth are often indication of misassemblies" |

- The pileup is vectorised in numpy: 5.4 Mb of read bases takes about 0.3 s per contig, and results are cached.

**Template display** (the **Template** tab above the editor, or `t`)
- Gap5's template plot, with the same colour scheme and with Y = template size on a log scale. On that scale, 500 bp Illumina inserts and 10–30 kb long reads separate naturally, much as Gap5's "Stacking Y Size" separates short- from long-read data.
- The shaded band is the library's normal insert range, estimated from the data.
- The depth track underneath shows read depth and template depth.
- Wheel to zoom, drag to pan, click to open the editor there.
- Each read's strand marker in the editor uses Gap5's names-panel colours for pair status.

**Read pairs as evidence**
- The library insert size comes from FR pairs (median ± 4 robust SD).
- Pairs whose mates sit in different contigs near the contig ends become join candidates (`pairs` evidence), with gap = insert − (distance of each read from its contig end).
- Outward-facing pairs across the two ends of one contig are counted as circular "wrap" evidence.

## How the evidence is scored (`evidence.py`)

- **graph link**: +3.
- **long-read bridges**: consecutive hits on one read leave contig X through one end and enter contig Y through another. Score +2·log2(1+n). Gap = read distance minus the unaligned contig ends. A read that leaves mid-contig is counted as a conflict, not as support.
- **read pairs**: mates in different contigs near the ends. Score +1.5·log2(1+n), with a gap estimate from the insert size.
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
- The clip consensus comes from iterated star alignment (edlib). On simulated ONT reads at ~6% error, 15–20 clips give a consensus that is more than 99% correct. Real ONT errors are more systematic (homopolymers), so a POA consensus would be more robust. The *Soft-clip consensus* box in the side panel still uses the older naive column vote.
- Gap filling only writes Ns. The aligned clip consensus is shown in the gap, but it is not yet spliced into the exported sequence.
- After each edit the working assembly is rebuilt. Unchanged contigs are reused, so this is quick for a bacterial genome (about 1–2 s) but not instant for large genomes.
- There is no base-level editing. Polishing belongs to dedicated tools.

## References

- Bonfield JK, Whitwham A (2010) Gap5—editing the billion fragment sequence assembly. *Bioinformatics* 26:1699–1703. doi:10.1093/bioinformatics/btq268
- Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094–3100. doi:10.1093/bioinformatics/bty191
- Wick RR et al. (2015) Bandage: interactive visualization of de novo genome assemblies. *Bioinformatics* 31:3350–3352. doi:10.1093/bioinformatics/btv383
- Cheng H et al. (2021) Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nat Methods* 18:170–175. doi:10.1038/s41592-020-01056-5
- Hunt M et al. (2015) Circlator: automated circularization of genome assemblies using long sequencing reads. *Genome Biol* 16:294. doi:10.1186/s13059-015-0849-0 (an established standalone tool for the same closure/rotation step)
- Kolmogorov M et al. (2019) Assembly of long, error-prone reads using repeat graphs (Flye). *Nat Biotechnol* 37:540–546. PMID 30936562
