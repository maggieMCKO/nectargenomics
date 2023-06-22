#!/usr/bin/env python3
# This script takes fasta file and checks if required sequence is in frame;
# if seq is not specified, checks each sequence;
# Outputs string provided by user

import argparse
import pyfastx


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', type=str, help='input fasta file')
parser.add_argument('-s', '--sequence', type=str, default='', help='name of teh sequence to check')
parser.add_argument('-n', '--name', type=str, help='the output string')
args = parser.parse_args()

# initiate list of lengths of sequnces in the fasta
seq_lens = []
fasta_dict = {}
for name, seq in pyfastx.Fasta(args.fasta, build_index=False):
    if args.sequence == '':
        seq_lens.append(len(seq))
    else:
        fasta_dict[name] = seq

if args.sequence == '':
    if all(x % 3 == 0 for x in seq_lens):
        print(args.name)
else:
    if len(fasta_dict[args.sequence]) % 3 == 0:
        print(args.name)
