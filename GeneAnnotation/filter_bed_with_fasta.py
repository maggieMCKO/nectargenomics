#!/usr/bin/env python
#

"""
This script takes bed/genePred file with predicted gene models
and a fasta file with filtered proteins;
it goes through the bed/genePred file and filters out entries not present in the fasta
"""

import argparse
from Bio import SeqIO
from collections import defaultdict

__author__ = "Ekaterina Osipova, 2019."


parser = argparse.ArgumentParser()
parser.add_argument('-a', '--anno', type=str, help='annotation file to filter in bed/gp format')
parser.add_argument('-f', '--fasta', type=str, help='fasta file with predicted gene models')
parser.add_argument('-gp', '--gp', action='store_true', help='if specified, expects annotation in genePred format')
args = parser.parse_args()


## Go through fasta and make a list of headers
fastaEntries = []
sequences = SeqIO.parse(open(args.fasta), 'fasta')
for entry in sequences:
    fastaEntries.append(entry.id)


## Go through the annotation file, filter out all entries that are not in the fastaEnries list
with open(args.anno, 'r') as inf:
    for line in inf.readlines():
        if args.gp:
            if line.split('\t')[0] in fastaEntries:
                print(line.rstrip())
        else:
            if line.split('\t')[3] in fastaEntries:
                print(line.rstrip())




