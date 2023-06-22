# Protein-coding gene annotation pipeline

## Tools and scripts
Tools and scripts used in the annotation pipeline

### Tools:
```
GenomeThreader
TAMA-GO (optional)
BUSCO
Augustus v3.3.2
genometools
EVidenceModeler
TOGA
bedtools
blastp
samtools
Stringtie
TransDecoder (optional)
```
### Scripts:
```
align2hints.pl (augustus)
join_mult_hints.pl (augustus)
augustus.pl (augustus)
partition_EVM_inputs.pl (EVM)
write_EVM_commands.pl (EVM)
EVM_to_GFF3.pl (EVM)
getSplitsOutsideGenes.pl (David)
partitionHints.pl (David)
evmFormatBraker.pl (David)
evmFormatAln.pl (David)
evmFormatIsoSeq.pl (David)
IncludeStopsToGenePredCDS.pl (David) (TAMA only)
keepOnlyCDSgenes.pl (David) (TAMA only)
filter_fasta_with_blast.py (Katya)
filter_bed_with_fasta.py (Katya)
getOverlappingTranscripts.py (Katya)
getUniqTranscripts.py (Katya)
add_utrs_from_stringtie.py (Katya)
```
```
ALL_SPECIES_LIST=birds_to_annotate.lst
```
## Building pairwise alignments to reference (chains)

## Making TOGA projections

## Processing transcriptome data
### Assemble transcriptome with Stringtie
```
stringtie -p 16 -l $db -o $db.strg.gtf $db.sorted.bam
```

## Alignment of proteins
### Align proteins and cDNA sequences with GenomeThreader

prepare and split protein file
```
md split_all_protein
md results_gth/
awk 'BEGIN {n_seq=0;} /^>/ {if (n_seq%100==0) {file=sprintf("prot_seq_%d.faa",n_seq);} print >> file; n_seq++; next;} {print >> file; }' < $all_protein.aa.fas
```

prepare bash_run_gth.sh script

prepare gth jobList and submit it
```
for i in $(ls split_all_protein/); do echo "bash_run_gth.sh $i"; done > jobList_gth
bash jobList_gth
```
process results
```
find results_gth/ -name "*.gff3" -type f -size -1k -delete
for i in $(ls results_gth/*.gff3); do echo $i; done > list.files.to.merge
```
create a jobList for sorting gff3 files:
```
for i in $(cat list.files.to.merge); do echo "gt gff3 -sort -tidy $i > ${i}.sorted"; done > jobList_gt_sort
para make gt_sort jobList_gt_sort

awk '{print $1".sorted"}' list.files.to.merge > file && mv file list.files.to.merge
gt merge $(cat list.files.to.merge) > $db.merged.gth.gff3
```


## De novo gene prediction with Augustus
### Augustus single

### Prepare hints
produce RNA hints for augustus
```
bam2hints --in=$db.sorted.bam --out=hints.$db.rna.gff
```
produce protein hints for augustus
```
align2hints.pl --in=$db.merged.gth.gff3 --out=/hints.$db.gth.gff --prg=gth
```

merge RNA and protein hints
```
cat hints.$db.rna.gff hints.$db.gth.gff > hints.$db.rna.prot.gff
```

join multiple hints
```
cat hints.$db.rna.prot.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl > hints.joined.$db.rna.prot.gff
```

Split genome and hint files

get split coordinates
```
getSplitsOutsideGenes.pl $pathToTOGAresults/$db/final.bed $db > splits.outsideToga.bed
```
split fasta
```
md scaffolds
faSplit byName $genomePath/gbdb-HL/$db/$db.fa scaffolds/
```

split hints
```
partitionHints.pl splits.outsideToga.bed hints.joined.$db.rna.prot.gff
```

Prepare and run augustus jobs
```
md augustus_results
write_aug_split_jobs.bash $db > jobList_augustus
para make augustus jobList_augustus
```

Process augustus results
```
for i in $(ls augustus_results/*.gtf); do echo $i; done > list.result.files.$db
joingenes -f list.result.files.$db -o $db.augustus.join.gtf
```

covert gtf -> gff3 (for each $db)
```
gt gtf_to_gff3 -tidy Augustus/$db.augustus.join.gtf | gt gff3 -tidy -sort > Augustus/$db.augustus.join.gff3
```


## Augustus cgp
### Augustus comparative mode


Prepare a multiple genome alignment with Multiz
Prepare a phylogeny of your species with branch lengths

Split maf file in 1Mb chunks
```
md mafSplit_outside_genes_ali/
mafSplit $splitPos.1MB.bed mafSplit_outside_genes/ multi.maf
```

delete files with no alignments
```
find mafSplit_outside_genes/ -name "*.maf" -type f -size 1k -delete
```
Prepare cgp hint files
```
md cgp
```

Prepare TOGA hints
```
md cgp/source_M_hints
```
for each $ref species, for each $db run
```
bedToGenePred $pathToTOGAresults/final.bed stdout | genePredToGtf file stdin stdout | grep -P '\t(CDS|start_codon|stop_codon)\t' | gtf2gff.pl --printIntron --includeStopInCDS --out=cgp/source_M_hints/$db.$ref.toga.gff
```
gff -> hints; join multiple hints
```
md cgp/hints_all/
```

for each $db run
```
for ref in $LIST_OF_REF; do grep -P '\t(CDS|intron)\t' cgp/source_M_hints/$db.$ref.toga.gff | cut -f1-8 | perl -pne 's/$/\tsrc=M/' | sort -n -k 4,4 | sort -s -n -k5,5 | sort -s -n -k3,3 |sort -s -k 1,1; done | join_mult_hints.pl > cgp/hints_all/${i}.toga.hints; done
```

add ref annotaions themselves
```
for ref in $LIST_OF_REF; do bedToGenePred $REF_ANNO_FILE stdout | genePredToGtf file stdin stdout | grep -P '\t(CDS|start_codon|stop_codon)\t' | gtf2gff.pl --printIntron --includeStopInCDS --out=cgp/hints_all/$ref.$ref.toga.gff; done

for ref in $LIST_OF_REF; do grep -P '\t(CDS|intron)\t' cgp/hints_all/$ref.$ref.toga.gff | cut -f1-8 | perl -pne 's/$/\tsrc=M/' | sort -n -k 4,4 | sort -s -n -k5,5 | sort -s -n -k3,3 |sort -s -k 1,1 | join_mult_hints.pl > cgp/hints_all/$ref.$ref.toga.hints; done
```

Prepare RNAseq hints
run bam -> hints for each $db
```
bam2hints --intronsonly --in=$db.sorted.bam --out=cgp/hints_all/$db.introns.hints
```
run bam2wig -> wig2hints for each $db
```
bam2wig $db.sorted.bam | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 ==prune=0.1 --src=W --type=ep --radius=4.5 --pri=3 --strand="." > cgp/hints_all/$db.ep.hints
```

combine RNAseq hints
```
cd cgp/hints_all
cat $db.introns.hints $db.ep.hints > $db.rna.hints
```

Combine all hints in one file per $db
```
cd ..
md hints_cgp/
```
for each $db run
```
cat hints_all/$db.rna.hints hints_all/$db.toga.hints | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl > hints_cgp/$db.all.cgp.hints
```
for each $ref
```
cp hints_all/$ref.toga.hints hints_cgp/$ref.all.cgp.hints
```
Make cgp.db

make hint table
```
for i in $ALL_SPECIES_LIST; do echo -e "$i\t$FULL_PATH_TO/cgp/hints_cgp/$i.all.cgp.hints"; done > hint.tbl
```
make genomes table
```
for i in $(cat $ALL_SPECIES_LIST); do echo -e "$i\t$genomePath/gbdb-HL/$i/$i.fa"; done > genomes.tbl
```
make sure that each genomic fasta exists

make cgp.db (NB: takes ~40 mins for 14 species)
```
make_cgp_db.sh genomes.tbl hint.tbl cgp.db
```
Prepare extrinsic config file (copy from Augustus and modify for your species)

Prepare and run cgp jobs

prepare cgp jobs. $model - Augustus model, e.g. human, chicken or custom
```
prepare_cgp_jobs_copyDb_short.sh $FULL_PATH_TO/mafSplit_outside_genes/ $FULL_PATH_TO/tree.nh $model $FULL_PATH_TO/cgp/genomes.tbl cgp.db $FULL_PATH_TO/cgp/$cgp_extrinsic.cfg $FULL_PATH_TO/cgp/output

for i in $(ls output/jobs/); do echo "sbatch output/jobs/$i"; done > jobList.sbatch.cgp
bash jobList.sbatch.cgp
```

Process results

make lists of files to merge
```
for db in $(cat $ALL_SPECIES_LIST); do echo $db;  find output/results/ -maxdepth 3 -name "$db.cgp.gff" > $db.list.to.sort; done
```
sort output gffs
```
for db in $(cat $ALL_SPECIES_LIST); do echo $db; for i in $(cat $db.list.to.sort); do gt gtf_to_gff3 -tidy $i | gt gff3 -tidy -sort -setsource cgp  > $i.sorted; done; done
```
merge sorted gffs
```
for db in $(cat $ALL_SPECIES_LIST); do awk -F "\t" '{$1=$1".sorted"; print}' $db.list.to.sort > $db.list.to.merge; done
for db in $(cat $ALL_SPECIES_LIST); do echo "gt merge \$(cat $db.list.to.merge) > $db.merged.cgp.gff3" ; done > jobList.merge.cgp
para make merge.cgp jobList.merge.cgp
```
Merge all de novo predictions
```
for db in $(cat $ALL_SPECIES_LIST); do gt merge Augustus/$db.augustus.join.gff3 cgp/$db.merged.cgp.gff3 > $db.all.denovo.gff3
```



## Filtering de novo predictions with blast
### Filter out repeats

on genome: get repeat files
```
cd /genome/gbdb-HL/$db
getRepeatClassOverview.perl $db -keepTempFiles
find . -maxdepth 1 -type f -regextype egrep -regex '.*[a-zA-Z0-9]{8}\.bed' | xargs -i mv {} repeats.$db.bed; done
find . -maxdepth 1 -type f -regextype egrep -regex '.*[a-zA-Z0-9]{8}\.bed.*' | g -v repeats | xargs -i rm {}
rm 2
```

back on falcon1: copy repeat files from genome
```
scp genome.pks.mpg.de:/genome/gbdb-HL/$db/repeats.* .
```

filter out repeats
```
bedtools intersect -v -f 0.1 -split -a $db.all.denovo.bed -b repeats.$db.bed -wa > ${db}_5.1/BlastFilter_repeats/$db.all.denovo.norepeats.bed
```
Run blast against swissprot database

bed -> gp -> prot
```
bedToGenePred $db.all.denovo.norepeats.bed stdout | genePredToProt -includeStop -starForInframeStops stdin $genomePath/gbdb-HL/$db/$db.2bit $db.all.denovo.norepeats.aa.fas
```
split prot.fastas
```
md split_fasta/
md results_blast/
cd split_fasta/
faSplit sequence ../$db.all.denovo.norepeats.aa.fas 300 split_fasta
```
prepare and run blastp jobList
```
for i in $( ls split_fasta/); do echo "blastp -evalue 1e-10 -num_threads 24 -db /projects/genome-bat/uniref/swissprot -query split_fasta/$i -outfmt 6 -out results_blast/${i%.fa}.BlastHits_out6.tsv"; done > jobList.batch.blastp
para make batch.blastp jobList.batch.blastp
```
merge blast results
```
cat results_blast/*.tsv > merged_${db}.denovo.norepeats.BlastHits_out6.tsv
```
Filter gene predictions with blast results

prepare species.dictionary.lst file with allowed species, like:
NOTVI Notophthalmus viridescens
HOPLI Hoplosternum littorale
DEMVE Demansia vestigiata
BUNCA Bungarus candidus
CHROW Chrotogale owstoni

filter for vertebrate hits OR >=200AA length
```
filter_fasta_with_blast.py -l 200 -b merged_${db}.denovo.norepeats.BlastHits_out6.tsv -f $db.all.denovo.norepeats.aa.fas -s vertebrates.dictionary.lst > $db.filtered.all.denovo.norepeats.aa.fas

filter_bed_with_fasta.py -a $db.all.denovo.norepeats.bed -f $db.filtered.all.denovo.norepeats.aa.fas > $db.filtered.all.denovo.norepeats.bed
```



## Generating consensus gene models
### Combine all evidence into consensus gene models with EVidenceModeler

Prepare input
format everything for EVM
```
md Evidences
```
format TOGA projections
run for each of $ref species
```
bedToGenePred $pathToTOGAresults/final.bed stdout | genePredToGtf file stdin stdout | gt gtf_to_gff3 -tidy | gt gff3 -tidy -sort -setsource $ref | evmFormatIsoSeq.pl > Evidences/$db.evm.toga.$ref.gff3
```
format gth
```
evmFormatAln.pl $db.merged.gth.gff3 > Evidences/$db.evm.protein.gff
```
format Stringtie
```
gt gtf_to_gff3 -tidy $db.strg.gtf | gt gff3 -tidy -sort -setsource stringtie | evmFormatIsoSeq.pl > Evidences/$db.evm.strg.gff3
```
format RNAseq processed with TAMA-GO (if you did this filtering step)
```
g 'full_length' $db.withCDS.strg.bed | g 'prot_ok' | perl -pne 's/\h+/\t/g' | perl -lane '@arr=split(/;/,$F[3]);$F[3]=$arr[0];print join("\t",@F)' \
| bedToGenePred stdin stdout | IncludeStopsToGenePredCDS.pl | genePredToGtf file stdin stdout | gt gtf_to_gff3 -tidy | gt gff3 -tidy -sort -setsource stringtie \
| evmFormatIsoSeq.pl | keepOnlyCDSgenes.pl > Evidences/$db.evm.strg.gff3
```
format cgp 
```
$db.merged.cgp.gff3 | gt gff3 -tidy -sort -setsource cgp | evmFormatBraker.pl | grep -v '' > Evidences/$db.evm.cgp.gff3
```
format augustus single
```
gt gtf_to_gff3 -tidy $db.augustus.join.gtf | gt gff3 -tidy -sort -setsource augustus | evmFormatBraker.pl | grep -v '' > Evidences/$db.evm.augustus.gff3
```
merge all genePreds for EVM
```
gt merge Evidences/$db.evm.toga.* Evidences/$db.evm.strg.gff3 Evidences/$db.evm.denovo.gff3 > Evidences/$db.evm.genepreds.gff3
```
run EVM
```
md EVM
cd EVM
```
prepare weights.txt file:
```
ABINITIO_PREDICTION denovo 1
OTHER_PREDICTION chicken 8
OTHER_PREDICTION greattit 8
OTHER_PREDICTION zebrafinch 8
OTHER_PREDICTION stringtie 2
PROTEIN gth 2
```
prepare fasta files for your $db
```
twoBitToFa $genomePath/gbdb-HL/$db/$db.2bit $genomePath/gbdb-HL/$db/$db.fa
```
prepare partitions
```
partition_EVM_inputs.pl --genome $genomePath/gbdb-HL/$db/$db.fa --gene_predictions Evidences/$db.evm.genepreds.gff3 --protein Evidences/$db.evm.protein.gff3 --segmentSize 1000000 --overlapSize 150000 --partition_listing evm_partitions_list.out
```
write EVM commands
```
write_EVM_commands.pl --genome $genomePath/gbdb-HL/$db/$db.fa --weights weights.txt --gene_predictions Evidences/$db.evm.genepreds.gff3 --protein Evidences/$db.evm.protein.gff3 --output_file_name $db.evm.out --partitions evm_partitions_list.out > evm_${db}_commands.jobList
```
split jobs (if $db has large number of scaffolds)
```
md batch_evm_${db}
shuf evm_${db}_commands.jobList | splitFile stdin 100 batch_evm_${db}/batch_
ls batch_evm_${db}/* | xargs -i echo "bash {}" >> batch.evm.$db
para make batch.evm batch.evm.$db
```
Process EVM results
covert results to gtf
```
for i in $(find . -name \*evm.out | awk -F "/" 'BEGIN { OFS="/";} {$NF=""; print $0 }'); do echo "EVM_to_GFF3.pl $i/$db.evm.out $i | gt gff3 -tidy -sort | gt gff3_to_gtf | perl formatEvmOutGtf.pl > $i/$db.evm.out.gtf"; done > jobList.evm2gtf

md batch_evm2gtf
shuf jobList.evm2gtf | splitFile stdin 100 batch_evm2gtf/batch_
ls batch_evm2gtf/* | xargs -i echo "bash {}" >> batch.evm2gtf
para make evm2gtf batch.evm2gtf

find . -name \*evm.out.gtf > gtf.evm.out.lst
joingenes -f gtf.evm.out.lst -o $db.evm.join.gtf
```
clean ./ and / after joingenes if needed
```
sed 's/\.\///g' $db.evm.join.gtf | sed 's/\///g' > file && mv file $db.evm.join.gtf
```
gtf -> gff3 -> gp -> bed
```
gt gtf_to_gff3 $db.evm.join.gtf | gt gff3 -tidy -sort -fixregionboundaries | gff3ToGenePred stdin $db.evm.join.gp
genePredToBed $db.evm.join.gp $db.evm.join.bed
```



#Adding high-confidence predictions back
TOGA projections shared between at least 2 reference species and TOGA projections of APPRIS principal

get uniq overlaps for each pair of reference TOGAs
```
$ref1
$ref2
getOverlappingTranscripts.py -f $PATH_TOGA_results_ref1/final.bed $PATH_TOGA_results_ref2/final.bed > uniq.overlap.$ref1.$ref2.bed
```
add ref TOGAs and APPRIS principal projections
```
getUniqTranscripts.py -f <(cat $db.evm.join.bed uniq.overlap* $PATH_TOGA_APPRIS/final.bed) > $db.added_togas.appris.bed 2> err
```


## Prediction of UTRs from transcriptome data
```
add_utrs_from_stringtie.py -a ${db}_5.1/$db.added_togas.appris.bed -r ${db}_5.1/RNAseq/$db.strg.bed > Add_UTRs/utrs.$db.added_togas.appris.bed
```


## Assessment of gene annotation completeness
assess annotation completeness with BUSCO

prepare template.busco.job
prepare protein (bed -> gp -> prot)
```
bedToGenePred $db.added_togas.appris.bed stdout | genePredToProt -includeStop -starForInframeStops stdin $genomePath/gbdb-HL/$db/$db.2bit $db.added_togas.appris.aa.fas
```
run busco
```
md busco
sbatch template.busco.job $db
```


## Assigning gene names
for each db in $ALL_SPECIES_LIST
```
bash_assign_gene_names.sh $db
```
For transcripts where gene name was not assigned, use blastp
get not_assigned transcripts
```
md Assign_from_blast/$db
grep -i 'rna\|ENSGAL' geneNames.$db.added_togas.appris.gp > Assign_from_blast/$db/not_assigned.$db.added_togas.appris.gp
cd Assign_from_blast/
```
prepare prot
```
genePredToProt -includeStop -starForInframeStops $db/not_assigned.$db.added_togas.appris.gp $genomePath/gbdb-HL/$db/$db.2bit $db/not_assigned.$db.added_togas.appris.aa.fas
```
split fasta
```
md $db/split_fasta
md $db/results_blast
cd $db/split_fasta
faSplit sequence ../not_assigned.$db.added_togas.appris.aa.fas 100 split_fasta
cd Assign_from_blast/
```
prepare and run blastp joblist
```
for i in $(ls $db/split_fasta/); do echo "blastp -evalue 1e-10  -num_threads 24 -db swissprot -query $db/split_fasta/$i -outfmt 6 -out $db/results_blast/${i%.fa}.BlastHits_out6.tsv"; done > jobList.batch.blastp
para make blastp jobList.batch.blastp
```
merge blast results
```
cat $db/results_blast/*.tsv > $db/merged_${db}.not_assigned.BlastHits_out6.tsv
```
find the best blast hit;
assign gene names to transcript names from vertebarate uniprot corresponding to the best blast hit
```
assign_genes_to_hits.py -u uniprot_sprot.fasta -b <(filter_blast_hits.py -b $db/merged_${db}.not_assigned.BlastHits_out6.tsv -s vertebrates.dictionary.lst -n 1 )  -s like > transc.gene.dict.$db.csv

cd ..
```
replace transcript names with uniprot gene names where possible
```
renameToHLscaffolds.py -a  geneNames.$db.added_togas.appris.gp  -d Assign_from_blast/transc.gene.dict.$db.csv > ext.geneNames.$db.added_togas.appris.gp
```

