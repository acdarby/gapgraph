"""gapgraph: a graph-aware, Gap5-style genome finishing editor.

Evidence sources combined on top of an assembly graph (GFA):
  * graph links (GFA L-lines from hifiasm / Flye / Verkko / SPAdes ...)
  * soft-clipped read ends (short reads, ONT, PacBio CLR/HiFi) at contig ends
  * split / multi-contig long-read alignments (PAF or SAM/BAM supplementary)
  * Hi-C contacts between contig ends (pairs format)
"""
__version__ = "0.1.0"

END_L = "L"  # 5' / head of a contig in its stored orientation
END_R = "R"  # 3' / tail
