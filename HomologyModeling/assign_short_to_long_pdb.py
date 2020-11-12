#!/usr/bin/env python3
#
import argparse
import pyfastx
import re
import sys


__author__ = "Ekaterina Osipova, 2020."


def read_fasta_dict(fasta_file):
    ## Reads alignment fasta into a dictionary

    fasta_dict = {}
    for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fasta_dict[name] = seq
    return fasta_dict


def find_long_sort(fasta_dict):
    ## Finds long and short seq

    long_seq = ''
    short_seq = ''
    for k, v in fasta_dict.items():
        if k.endswith('long'):
            long_seq = v
        elif k.endswith('short'):
            short_seq = v
    # check if we got something
    if long_seq == '' or short_seq == '':
        print('ERROR:Did not find long or short sequence!')
        sys.exit(1)
    return long_seq, short_seq


def compare_two_seq(seq_long, seq_short):
    ## Makes correspondance between indexes in long and short sequences

    i = 1
    j = 1
    index_to_keep = {}
    for k in range(len(seq_long)):
        if seq_long[k] != '-':
            i += 1
        if seq_short[k] != '-':
            j += 1
        if (seq_long[k] != '-') and (seq_short[k] != '-'):
            index_to_keep[i] = j
    return index_to_keep


def read_replacement(file):
    ## Read replacement file into a list

    replace_values = []
    with open(file, 'r') as inf:
        for line in inf.readlines():
            replace_values.append(float(line.rstrip()[:4]))
    return  replace_values


def replace_pdb_bfactor(pdb, replace_values, index_to_keep):
    ## Reads PDB file, replace b-factor column with values from the list

    with open(pdb, 'r') as inf:
        for line in inf.readlines():
            if line.startswith('ATOM'):
                res_number = int(line.rstrip()[22:26].rstrip())
                old_pdb_line = line.rstrip()[:60]
                last_field = line.rstrip()[-12:]
                if res_number in index_to_keep:
                    new_index = index_to_keep[res_number]
                    new_bfactor = replace_values[new_index - 1]
                else:
                    new_bfactor = 0.00
                # output new pdb line
                print('{} {:.2f}{}'.format(old_pdb_line, new_bfactor, last_field))


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, help='fasta file with protein alignment')
    parser.add_argument('-p', '--pdb', type=str, help='file in pdb format with only ATOM entries')
    parser.add_argument('-r', '--replace', type=str, help='file with a list of values to replace 11th column in pdb')
    args = parser.parse_args()

    ## Read fasta
    fasta_dict = read_fasta_dict(args.fasta)

    ## Find long and short seq
    long_seq, short_seq = find_long_sort(fasta_dict)

    ## Make a list of indexes of a long sequence to assign p-values to
    index_to_keep = compare_two_seq(long_seq, short_seq)

    ## Get PDB file and assign p-values to corresponding positions in the structure
    replace_values = read_replacement(args.replace)

    ## Read PDB file, replace 11th column with values from the list
    replace_pdb_bfactor(args.pdb, replace_values, index_to_keep)


if __name__ == "__main__":
        main()