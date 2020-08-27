#!/usr/bin/env python
# This script takes fasta file and checks if each sequence in it is multiple of three;
# Outputs string you provides it as argument

import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', type=str, help='input fasta file')
parser.add_argument('-n', '--name', type=str, help='the output string')
args = parser.parse_args()

# initiate list of lengths of sequnces in the fasta
seqLens = []
sequences = SeqIO.parse(open(args.fasta), 'fasta')
for entry in sequences:
    seqName, seq = entry.id, str(entry.seq)
    seqLens.append(len(seq))

#print(seqLens)
if all(x % 3 == 0 for x in seqLens):
        print(args.name)
