Current data show no signal of Ebola virus adapting to humans
=============================================================
Stephanie J. Spielman, Austin G. Meyer, Claus O. Wilke

This repository contains all data and code to reproduce the analysis.



## Naming convention for alignments

* `.aln`: alignment, fasta format
* `_unique.aln`: alignment, fasta format, only uniques retained
* `_w_outgroup.aln`: alignment, fasta format, with appropriate outgroup
* `_w_outgroup_unique.aln`: as `_w_outgroup.aln` but only uniques retained

Note that we need to distinguish between duplicate strains (as in two sequences from the same patient) and duplicate sequences (two sequences share no mutations). We will use the word "duplicates" or "no-duplicates" to refer to duplicate strains, and the word "uniques" to refer to sequences with at least one difference.
