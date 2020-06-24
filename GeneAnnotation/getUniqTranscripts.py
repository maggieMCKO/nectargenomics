#!/usr/bin/env python
#
#
# To get unique transcripts of a file run e.g:\
# e.g usage: getUniqTrascripts.py -f HLparMaj1.ncbi.bed12 > uniq.HLparMaj1.ncbi.bed12 2> dupl.HLparMaj1.ncbi.bed12
#

import argparse
from collections import defaultdict
import sys


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filebed', type=str, help='bed12 file')
parser.add_argument('-c', '--color', action='store_true', help='if specified, RGB column must also be identical in duplicated transcripts')
parser.add_argument('-s', '--score', action='store_true', help='if specified, score column must also be identical in duplicated transcripts')
args = parser.parse_args()

# make initial dictionary of unique transcripts
overlapingTranscripts = defaultdict(list)
with open(args.filebed, 'r') as inf:
    for line in inf.readlines():
        transcrName, score, color = line.split()[3], line.split()[4], line.split()[8]

        # get transcript info considering which fields are important
        if args.color and args.score:
            # both score and color are considered
            transcrInfo = '\t'.join(line.split()[:3] + line.split()[4:])
        elif args.color and not args.score:
            # only color is considered
            transcrInfo = '\t'.join(line.split()[:3] + line.split()[5:])
        elif args.score and not args.color:
            # only score is considered
            transcrInfo = '\t'.join(line.split()[:3] + line.split()[4:8] + line.split()[9:])
        else:
            # both color and score are ignored
            transcrInfo = '\t'.join(line.split()[:3] + line.split()[5:8] + line.split()[9:])


        if (transcrInfo in overlapingTranscripts):
            sys.stderr.write('DUPLICATION: {}\t{}'.format(transcrName, '\t'.join(transcrInfo.split())+'\n'))
        else:
            overlapingTranscripts[transcrInfo].append((transcrName, score, color))

# output unique elements of the overlapingTranscripts dictionary
for transcrInfo in overlapingTranscripts:
    transcrNames = overlapingTranscripts[transcrInfo]

    # assemble a new annotation line considering fields excluded before
    if args.color and args.score:
        # both score and color are considered
        newBedLine = '\t'.join(transcrInfo.split()[:3]) + '\t' + transcrNames[0][0] + '\t' + \
                     transcrNames[0][1] + '\t' + '\t'.join(transcrInfo.split()[4:7]) + '\t' + \
                     transcrNames[0][2] + '\t' + '\t'.join(transcrInfo.split()[8:])
    elif args.color and not args.score:
        # only color is considered
        newBedLine = '\t'.join(transcrInfo.split()[:3]) + '\t' + transcrNames[0][0] + '\t' + \
                     transcrNames[0][1] + '\t' + '\t'.join(transcrInfo.split()[3:6]) + '\t' + \
                     transcrNames[0][2] + '\t' + '\t'.join(transcrInfo.split()[7:])
    elif args.score and not args.color:
        # only score is considered
        newBedLine = '\t'.join(transcrInfo.split()[:3]) + '\t' + transcrNames[0][0] + '\t' + \
                     transcrNames[0][1] + '\t' + '\t'.join(transcrInfo.split()[4:7]) + '\t' + \
                     transcrNames[0][2] + '\t' + '\t'.join(transcrInfo.split()[7:])
    else:
        # both color and score are ignored
        newBedLine = '\t'.join(transcrInfo.split()[:3]) + '\t' + transcrNames[0][0] + '\t' + \
                     transcrNames[0][1] + '\t' + '\t'.join(transcrInfo.split()[3:6]) + '\t' + \
                     transcrNames[0][2] + '\t' + '\t'.join(transcrInfo.split()[6:])


    print(newBedLine)
