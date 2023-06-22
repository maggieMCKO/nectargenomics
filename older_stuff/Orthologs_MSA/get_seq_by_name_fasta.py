#!/usr/bin/env python3
# This script takes fasta file and extracts a specified fasta sequence by name;
# Output to stdout

import argparse
import pyfastx


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', type=str, help='input fasta file')
parser.add_argument('-n', '--names', type=str, help='names of the sequences (comma-separated) to extract')
parser.add_argument('-p', '--prefix', type=str, default='', help='prefix you want to give your fasta name')
parser.add_argument('-s', '--suffix', type=str, default='', help='suffix you want to give your fasta name')
parser.add_argument('-ab', '--allbut', action='store_true', help='if specified, outputs every sequnce BUT specified in names')
args = parser.parse_args()


names = args.names.split(',')
for seq_name, seq in pyfastx.Fasta(args.fasta, build_index=False):    
    
    if args.allbut:
        if seq_name not in names:
            print(">" + args.prefix + seq_name + args.suffix)
            print(seq)
    else:       
        if seq_name in names:
             print(">" + args.prefix + seq_name + args.suffix)
             print(seq)
