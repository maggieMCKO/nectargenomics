#!/usr/bin/env python
#

"""
This script takes files in bed12 format (at least TWO!) and outputs overlapping unique transcripts;
it assigns transcript names using the first specified file.
Duplicated transcripts are going to stderr

To get overlaping transcripts of several annotations run e.g: \
getOverlapingTrascripts.py -f uniq.HLfloFus1.galGal6.ncbi.bed uniq.HLfloFus1.HLtaeGut4.ncbi.bed > uniq.overlap.HLfloFus1.galGal6.HLtaeGut4.bed
"""

import argparse
from collections import defaultdict
import sys

__author__ = "Ekaterina Osipova, 2019."



def read_all_annotations(bedList):
    ## Reads all transcripts of all annotation files into 2D list

    # initiate a 2D list to store lists of elements from each file
    allTranscripts_allFiles = []

    # initiate a dictionary to store the part of a transcript that is allowed to vary
    transcripts_dict = defaultdict(list)

    for bed in bedList:
        with open(bed, 'r') as inf:
            transcripts_oneFile = []
            for line in inf.readlines():
                transcrName, score, color = line.split()[3], line.split()[4], line.split()[8]
                transcrInfo = '\t'.join(line.split()[:3] + line.split()[5:8] + line.split()[9:])

                # store essential info for each transcript of this file in the list
                transcripts_oneFile.append(transcrInfo)

                # at the same time keep corresponding variable part of a transcript in the dictionary
                transcripts_dict[transcrInfo].append((transcrName, score, color))

            # add essential (that suppose to be identical) transcript info to the allTranscrpts 2D list
            allTranscripts_allFiles.append(transcripts_oneFile)
    return allTranscripts_allFiles, transcripts_dict


def print_overlapping_transcripts(overlap_transcripts, transcripts_dict):
    ## Outputs elements of the overlap_transcripts getting first available variable part from the transcripts_dict

    for transcrInfo in overlap_transcripts:
        transcrVars = transcripts_dict[transcrInfo]
        # combine a new line
        newBedLine = '\t'.join(transcrInfo.split()[:3]) + '\t' + transcrVars[0][0] + '\t' + \
                     transcrVars[0][1] + '\t' + '\t'.join(transcrInfo.split()[3:6]) + '\t' + \
                     transcrVars[0][2] + '\t' + '\t'.join(transcrInfo.split()[6:])
        print(newBedLine)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filelist', nargs='*', type=str, help='list of bed12 files to extract identical\
                                                          overlapping transcripts from. Needs at least TWO files')
    args = parser.parse_args()

    ## Read all annotation files into 2D list
    allTranscripts_allFiles, transcripts_dict = read_all_annotations(args.filelist)

    ## Get overlapping transcripts
    overlap_transcripts = set(allTranscripts_allFiles[0]).intersection(*allTranscripts_allFiles)

    ## Output overlapping transcripts
    print_overlapping_transcripts(overlap_transcripts, transcripts_dict)


if __name__ == '__main__':
    main()

