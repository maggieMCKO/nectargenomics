### Reference species annotation
- GCF_000002315.6_GRCg6a_genomic.gff
- To extract 4d codons, we need to have matching chromosome names in the alignment (.maf) and annotation (.gff), and the alignment has to be on the same chromosome.
- Need to split both annotation and alignment by chromosome.

#### 1. Prepare .gff
1. Replace chromosome names into this format: 'galGal6.chr1'
2. Split up the gff by chromosome.

#### 2. Prepare .maf
1. Use 'maffilter' to split by chromosome.
2. Run 'maffilter_selchr_SUBMIT_de.sh' which calls and modifies maffilter_selchr_slurm_de.sh and optionfile_tmp.maffilter according to each chromosome
