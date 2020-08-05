#!/bin/bash
#SBATCH -J gth.HLmalCya1               # Job name
#SBATCH -n 1                 # Run a single task        
#SBATCH -c 1                 # Number of CPU cores per task
#SBATCH --time=8:00:00	     # Set walltime limit
#SBATCH --mem-per-cpu=20000  # Job memory request
#SBATCH -o log.gth.HLmalCya1     # Standard output and error log

## report what node we're on
echo "Working on $HOSTNAME" 

## make working dir in /tmp
mkdir -p /tmp/$USER/$SLURM_JOB_ID/
echo "Building in /tmp/$USER/$SLURM_JOB_ID/"

## goto working dir
cd /tmp/$USER/$SLURM_JOB_ID/

## put genome here 
cp /projects/hillerlab/genome/gbdb-HL/HLmalCya1/HLmalCya1.fa .

## put proteins here
cp /projects/project-osipova/NectarivoryProject/GenomeAnnotation/Old_annotation_runs/HLfloFus1_1.0/Protein_evidence/split_all_protein/$1 .

## gth command
gth -skipalignmentout -gff3out -paralogs -prseedlength 20 -prminmatchlen 20 -prhdist 2 -minmatchlen 32 -seedlength 32 -gcmincoverage 70 -species chicken -genomic HLmalCya1.fa -protein $1 -o HLmalCya1.$1.gthAln.gff3

## move result file back to project space

mv HLmalCya1.$1.gthAln.gff3 /projects/project-osipova/NectarivoryProject/GenomeAnnotation/HLmalCya1_5.1/Protein_evidence/results_gth/

## clean up on the node! no matter what..
rm -r /tmp/$USER/$SLURM_JOB_ID/
