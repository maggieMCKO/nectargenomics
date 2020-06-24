#!/bin/bash

db=$1
model="chicken"
wdir=/projects/project-osipova/NectarivoryProject/GenomeAnnotation/${db}_5.1/
#fhsplit=$wdir/Augustus/hints.joined.$db.rna.prot.gff
splits=$wdir/splits.outsideToga.bed
hsplit=$wdir/Augustus/hints_split/
outdir=$wdir/Augustus/augustus_results


while read split;do
	chr=$(echo $split | awk '{print $1}')
	start=$(echo $split | awk '{print $2}')
	end=$(echo $split | awk '{print $3}')
	if [ -f $hsplit/$chr.$start.$end.split.hints ]; then
	
		echo "augustus --species=$model --extrinsicCfgFile=/projects/genome-bat/.batcave/mybin/augustus-3.3.1/config/extrinsic/extrinsic.M.RM.E.W.P.PB.cfg --allow_hinted_splicesites=gcag,atac --codingseq=on --hintsfile=$hsplit/$chr.$start.$end.split.hints --AUGUSTUS_CONFIG_PATH=/projects/genome-bat/.batcave/mybin/augustus-3.3.1/config/ --alternatives-from-evidence=true --predictionStart=$start --predictionEnd=$end --UTR=off $wdir/Augustus/scaffolds/$chr.fa > $outdir/$chr.$start.$end.pred.gtf"
	else
		echo "augustus --species=$model --extrinsicCfgFile=/projects/genome-bat/.batcave/mybin/augustus-3.3.1/config/extrinsic/extrinsic.M.RM.E.W.P.PB.cfg --allow_hinted_splicesites=gcag,atac --codingseq=on --AUGUSTUS_CONFIG_PATH=/projects/genome-bat/.batcave/mybin/augustus-3.3.1/config/ --alternatives-from-evidence=true --predictionStart=$start --predictionEnd=$end --UTR=off $wdir/Augustus/scaffolds/$chr.fa > $outdir/$chr.$start.$end.pred.gtf"
	fi
done < $splits
