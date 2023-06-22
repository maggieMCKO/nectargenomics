#!/bin/bash

set -eo pipefail

if [ $# -lt 3 ];then
	echo "Not enough arguments! See usage"
	echo ""
	echo "USAGE: $0 [Genome Table File] [Hint Table File] [Output DB name] [Working directory, cannot be on Lustre!]"
	exit 1
fi

## workdir should be in /tmp/
WORKDIR=$4
genomeTbl=$1
hintTbl=$2
name=$3

if [ "$WORKDIR" = "" ];then
	WORKDIR=/tmp/$USER/
fi

if [ ! -d $WORKDIR ];then
	mkdir $WORKDIR
fi

echo "Building db $name in $WORKDIR"

while read genLine;do 
spec=$(echo "$genLine" | cut -f1)
echo $spec
genFile=$(echo "$genLine" | cut -f2)
echo "---------------"
echo $geneFile
echo "Running command: /projects/genome-bat/.batcave/mybin/Augustus/bin/load2sqlitedb --noIdx --species=$spec --dbaccess=$WORKDIR/$name $genFile"
/projects/genome-bat/.batcave/mybin/Augustus/bin/load2sqlitedb --noIdx --species=$spec --dbaccess=$WORKDIR/$name $genFile; 
done < $genomeTbl

while read hintLine;do 
spec=$(echo "$hintLine" | cut -f1)
hintFile=$(echo "$hintLine" | cut -f2)
/projects/genome-bat/.batcave/mybin/Augustus/bin/load2sqlitedb --noIdx --species=$spec --dbaccess=$WORKDIR/$name $hintFile; 
done < $hintTbl

/projects/genome-bat/.batcave/mybin/Augustus/bin/load2sqlitedb --makeIdx --dbaccess=$WORKDIR/$name

sqlite3 -header -column $WORKDIR/$name "SELECT count(*) AS '#hints',typename,speciesname FROM (hints as H join featuretypes as F on H.type=F.typeid) natural join speciesnames GROUP BY speciesid,typename;"

mv $WORKDIR/cgp.db .
