#!/bin/bash
set -euo pipefail


#################################################
#					            				#
# Script to prepare Augustus CGP jobs and files #
#                                               #
#################################################

if [ $# -lt 7 ];then
	echo "Not enough arguments! See usage"
	echo ""
	echo "USAGE: $0 [Directory with split maf files] [Phylogenetic tree with neutral branch lengths] [Species config to use] [Genome Table File] [Database with hints and genomes] [Extrinsic config file] [Output directory]"
	exit 1
fi

## directory containing split maf files
mafDir=$1

## tree with neutral branch lengths for maf species
TREE=$2

## species to use for Augustus. Mammals -> human; Birds -> chicken;
model=$3

## genome table file as per cgp manual
genomeTbl=$4

## database with hints and genome offsets. Create with make_cgp_db.sh script 
db=$5

## extrinsic config file
config=$6

## name of the output dir for all files
outDir=$7

## get database file name
dbName=`basename $db`

## make output directories
mkdir -p $outDir/jobs/
mkdir -p $outDir/logs/
mkdir -p $outDir/results/

## loop over mafs in split dir and write a batch script for each
for chunk in $(ls $mafDir/*maf); do
id=`basename $chunk .maf`
echo "CREATING JOB FOR $id"

jobScript=$(cat <<EOF
#!/bin/bash
#SBATCH -J cgp-$id       # Job name
#SBATCH -n 1                # Run a single task        
#SBATCH -c 1                # Number of CPU cores per task
#SBATCH -t 24:00:00        # Walltime limit
#SBATCH --mem-per-cpu=20000  # Job memory request
#SBATCH -o $outDir/logs/cgp-$id.log   # Standard output and error log

echo "Running on \$HOSTNAME"
WORKDIR=/tmp/$USER/\$SLURM_JOBID

echo "Building sqlite db \$WORKDIR/$dbName"
echo "Job commencing at \`date\`"
SECONDS=0
mkdir -p \$WORKDIR
mkdir -p $outDir/results/$id

echo "....."
echo "Copying sqlite db for genomes and hints to \$WORKDIR"
echo -e ".....\n\n"
cp $db \$WORKDIR

### Loaded genomes and hints into database

/projects/genome-bat/.batcave/mybin/Augustus/bin/augustus --species=$model --dbaccess=\$WORKDIR/$dbName --extrinsicCfgFile=$config --dbhints=true --treefile=$TREE --alnfile=$chunk --speciesfilenames=$genomeTbl --softmasking=1 --UTR=off --/CompPred/outdir=$outDir/results/$id/
augtest=\$?

if [ "\$augtest" -eq 0 ];then
	echo "Augustus-cgp completed successfully for $chunk"
	echo "Completed at \`date\`"
	echo "Took \$SECONDS"
else
	echo "Augustus-cgp FAILED for $chunk"
	echo "augustus --species=$model --dbaccess=\$WORKDIR/$dbName --extrinsicCfgFile=$config --dbhints=true --treefile=$TREE --alnfile=$chunk --speciesfilenames=$genomeTbl --softmasking=1 --UTR=off --/CompPred/outdir=$outDir/results/$id/" >> $outDir/FAILEDCOMMANDS
fi

## remove workdir ANYWAY!
rm -rf \$WORKDIR

EOF
)

echo "$jobScript" > $outDir/jobs/$id.job
done
