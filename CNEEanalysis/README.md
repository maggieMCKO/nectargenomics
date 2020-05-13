## phyloFit
- Reference species annotation: GCF_000002315.6_GRCg6a_genomic.gff
- To extract 4d codons, we need to have matching chromosome names in the alignment (.maf) and annotation (.gff), and the alignment has to be on the same chromosome.
- Need to split both annotation and alignment by chromosome.

### File preprocessing
#### 1. Prepare .gff
1. use script ```00_maf_gff_preprocessing/1_gff/1_fix_galgalanno_maffilter_includeGalgal.sh```
    - Replace chromosome names into this format: 'galGal6.chr1'
    - fix six rows with 'Curated  Genomic' which messed up the gff (one additional tab), replace it to 'Curated_Genomic'
    - Split up the gff by chromosome.

#### 2. Prepare .maf
1. Use script ```00_maf_gff_preprocessing/2_maf/maffilter_selchr_SUBMIT_de.sh``` which calls and modifies maffilter_selchr_slurm_de.sh and optionfile_tmp.maffilter according to each chromosome
    - Use 'selchr' function in 'maffilter' to split by chromosome.

### Running phyloFit
1. Use script ```01_phylofit/run_phylofit_perChr_de.sh```
2. Output: ```01_phylofit/nonconserved_4d.mod```

---
## phyloP
### File preprocessing
#### 1. Prepare feature.bed (in this case, the CNEEs)
1. Use script ```02_phylop/01_getcnee.sh```
#### 2. Running phyloFit to get the neutral model
- See above

### Running phyloP
1. Get an estimate of departures from neutrality for the CNEE elements
    - use script ```02_phylop/02_run_phylop_de_v0.sh```
2. Use --subtree or --branch functions
((under construction))
    - script ```02_phylop/02_run_phylop_subtree_de_v0.sh```


---
## phyloacc
### Preparation before running
#### 1. Convert the multiple alignment from .maf to .fasta format
- use script ```3_maf_to_fasta/maf_to_fasta_de_v0.sh```

#### 2. Running phyloFit to get a rooted tree with branch length
- use script ```../01_phylofit/run_phylofit_perChr_de.sh```
- output: ```../01_phylofit/nonconserved_4d.mod```

#### 3. CNEE coordinates (.bed)
- use script ```../02_phylop/01_getcnee.sh```

#### 4. Parameter file (species names and parameters for MCMC)
