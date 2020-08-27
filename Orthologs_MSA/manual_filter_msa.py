#!/usr/bin/env python3
#
# This script takes multiple alignment in fasta format and performs 2-step filtering:
# 1) moves codon by codon and checks how much seq align;
# 2) moves species by species and checks if enough seq present

import argparse
import pyfastx
import re
import numpy as np
from operator import itemgetter


__author__ = "Ekaterina Osipova, 2020."



def read_fasta(fasta_file, ref):
    ## Reads alignment fasta (pyfastx module is faster than biopython) into fasta_dict; checks correctness

    fasta_dict = {}
    for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fasta_dict[name] = seq

    # check if all sequences have the same length
    seq_lens = set([len(fasta_dict[i]) for i in fasta_dict])
    try:
        assert len(seq_lens) == 1
    except AssertionError as e:
            e.args += ('ERROR: sequences in fasta are not the same length!', 1)
            raise

    # # check if sequence len is multiple of 3
    # ref_name = [name for name in fasta_dict if ref in name][0]
    # try:
    #     assert ref_name != ''
    # except AssertionError as e:
    #     e.args += ('ERROR: could not find ref in fasta! {}'.format(ref), 1)
    #     raise
    # ref_seq = fasta_dict[ref_name]
    #
    # if len(ref_seq) % 3 != 0:
    #     fasta_dict = fix_fasta_by_ref(fasta_dict, ref_name)

    return fasta_dict


def fix_fasta_by_ref(fasta_dict, ref_name):
    ## Tries to fix MSAs not multiple of 3 by trimming the beginning till the first ref ATG

    START = 'ATG'
    start_index = fasta_dict[ref_name].upper().find(START)

    # trim all sequences in the dictionary
    fixed_fasta_dict = {}
    for name in fasta_dict:
        fixed_fasta_dict[name] = fasta_dict[name][start_index: ]
    return fixed_fasta_dict


def codons_to_bool(seq, good_codon_pattern):
    ## Converts sequence into a list of True/False: T - when codon is good, F - otherwise

    codons = [seq[i: i + 3] for i in range(0, len(seq), 3)]
    codon_bool_seq = [bool(re.match(good_codon_pattern, codon)) for codon in codons]
    return codon_bool_seq


def check_bool_matrix(matrix, min_good):
    ## Checks 2D matrix of bool values if each column has enough True(=good) values

    row_number = matrix.shape[0]
    column_number = matrix.shape[-1]
    positions_to_keep = []

    for i in range(column_number):
        column = list(matrix[:, i])
        if column.count(True) >= round(min_good * row_number):
            positions_to_keep.append(i)
    return positions_to_keep


def filter_site_by_site(fasta_dict, positions_to_keep, mask, good_codon_pattern, STOPS):
    ## Site by site filtering

    filtered_fasta_dict = {}
    for name in fasta_dict:
        seq = fasta_dict[name]
        codons = [seq[i: i + 3] for i in range(0, len(seq), 3)]
        if len(positions_to_keep) > 0:
            filtered_codons = list(itemgetter(*positions_to_keep)(codons))
        else:
            filtered_codons = []

        # mask bad codons and stop-codons with NNN
        if mask:
            filtered_codons = [codon if (bool(re.match(good_codon_pattern, codon)) and (codon not in STOPS)) else 'NNN' for codon in filtered_codons]

        filtered_fasta_dict[name] = ''.join(filtered_codons)
    return filtered_fasta_dict


def filter_seq_by_seq(fasta_dict, min_seqfraction, min_seqlen, good_codon_pattern):
    ## Sequence by sequence filtering

    filtered_names = []
    for name in fasta_dict:
        bool_codons = codons_to_bool(fasta_dict[name], good_codon_pattern)
        if (bool_codons.count(True) >= round(min_seqfraction * len(bool_codons))) \
                and (len(bool_codons) >= min_seqlen):
            filtered_names.append(name)
    return filtered_names


def print_alignment(fasta_dict, good_sequences):
    ## Outputs alignment to stdout

    for name in fasta_dict:
        if name in good_sequences:
            print('>{}'.format(name))
            print(fasta_dict[name])



def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--alignment', type=str, help='alignment in fasta format')
    parser.add_argument('-mc', '--mincodon', type=float, default=0.8, help='min codon fraction to align for a column to stay; default: 80%')
    parser.add_argument('-ms', '--minseq', type=float, default=0.3, help='min sequence fraction to align for it to stay; default: 30%')
    parser.add_argument('-ml', '--minaalen', type=int, default=50, help='min length of a sequence in AA to stay; default: 50')
    parser.add_argument('-m', '--mask', action='store_true', help='mask ambiguous codons and stop-codons with NNN')
    parser.add_argument('-r', '--ref', type=str, default='galGal6', help='ref species code to look for in fasta, default: galGal6')
    args = parser.parse_args()

    ## Define what is considered to be a good codon
    good_codon_pattern = r'[ATGCatgc]{3}'
    STOPS = ['TAG', 'TGA', 'TAA']

    ## Read fasta alignment into a dictionary
    fasta_dict = read_fasta(args.alignment, args.ref)

    ## Prepare bool matrix for the whole alignment
    bool_matrix_alignment = np.array([codons_to_bool(fasta_dict[fa], good_codon_pattern) for fa in fasta_dict])

    ## Identify good positions to keep
    positions_to_keep = check_bool_matrix(bool_matrix_alignment, args.mincodon)

    ## Filter each sequence in the alignment codon by codon using positions_to_keep list
    column_filtered_fasta_dict = filter_site_by_site(fasta_dict, positions_to_keep, args.mask, good_codon_pattern, STOPS)

    ## Filter out sequences with less than $fraction columns
    filtered_fasta_names = filter_seq_by_seq(column_filtered_fasta_dict, args.minseq, args.minaalen, good_codon_pattern)

    ## Output new alignment
    print_alignment(column_filtered_fasta_dict, filtered_fasta_names)


if __name__ == "__main__":
    main()