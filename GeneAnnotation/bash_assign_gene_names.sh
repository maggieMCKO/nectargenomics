#!/usr/bin/env bash
#


db=$1
ens_toga=/projects/project-osipova/NectarivoryProject/TOGA_galGal6_ensembl_iso/ensGene.galGal6.all.bed12
anno_bed=/projects/project-osipova/NectarivoryProject/TOGA_ref_species_ncbi/TOGA_galGal6_new_run/galGal6.ncbi.anno.bed12
out_file=geneNames.$db.added_togas.appris.gp
scripts_dir=/projects/project-osipova/PythonScripts/
echo $db

overlapSelect -statsOutput -strand -inCds -inFmt=bed -selectFmt=bed $ens_toga $anno_bed stdout | cut -f1-2 | rev | cut -d"." -f2- | rev | awk '{print $2"\t"$1}' > $db.ensTrans.ncbi.dictionary.txt

$scripts_dir/renameToHLscaffolds.py  -a $db.ensTrans.ncbi.dictionary.txt -d /projects/project-osipova/NectarivoryProject/GenomeAnnotation/ensTrans.geneName.galGal6.txt | awk '{print $2"\t"$1}' | sed 's/\t/,/g' > $db.geneName.ncbi.dictionary.txt

$scripts_dir/renameToHLscaffolds.py -a <(bedToGenePred $anno_bed stdout) -d $db.geneName.ncbi.dictionary.txt > $out_file
