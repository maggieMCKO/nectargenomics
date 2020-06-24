#!/usr/bin/env python3
#

"""
This script parses blast hits file (outfmt 6), filters hits of only allowed species, outputs N best hits
"""

import sys
import argparse
from collections import defaultdict

__author__ = "Ekaterina Osipova, 2020."



def read_species_file(species_file):
    ## Reads species file and makes a list of allowed species codes: [CHICK, URILO, ..]

    species_codes = []
    with open(species_file, 'r') as species_inf:
        for line in species_inf.readlines():
            species_codes.append(line.split('\t')[0])
    return species_codes


def filter_blast_hits(blast_hits, species_codes):
    ## Reads blast outfmt6 file into a dict; filters for allowed species codes
    ## qseqid, rseqid, pid, alilen, mism, gapop, qst, qend, rst, rend, eval, bitscore

    blast_hits_dict = defaultdict(list)
    with open(blast_hits, 'r') as inf:
        for line in inf.readlines():
            qseqid = line.split('\t')[0]
            ref_code = (line.split('\t')[1]).split('_')[-1]
            if ref_code in species_codes:
                entry_info = '\t'.join(line.split('\t')[1:]).rstrip()
                blast_hits_dict[qseqid].append(entry_info)
    return blast_hits_dict


def get_best_hits(blast_hits_dict, n):
    ## Gets best-scoring (considers bit-score) hits for each query entry

    best_hits_dict = {}
    for qseqid in blast_hits_dict:
        hits_list = blast_hits_dict[qseqid]
        hits_list.sort(key = lambda x: float(x.split()[-1]), reverse = True)
        if (n == 0) or (len(hits_list) < n):
            best_hits_dict[qseqid] = hits_list
        else:
            best_hits_dict[qseqid] = hits_list[ :n]
    return best_hits_dict


def print_best_hits(best_hits_dict):
    ## Outputs n best hits to stdout

    for qseqid in best_hits_dict:
        for hit in best_hits_dict[qseqid]:
            print('{}\t{}'.format(qseqid, hit))
    return


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--blasthits', type=str, help='blast hits output in format 6')
    parser.add_argument('-s', '--species', type=str, help='species file, latin names')
    parser.add_argument('-n', '--nbesthits', type=int, default=0, help='number of best hits to return; default: all (=0)')
    args = parser.parse_args()


    ## Read file with allowed species codes
    species_codes = read_species_file(args.species)

    ## Filter blast hits only for allowed species
    blast_hits_dict = filter_blast_hits(args.blasthits, species_codes)

    ## Get number of best hits specified by user
    best_hits_dict = get_best_hits(blast_hits_dict, args.nbesthits)

    ## Output best hits for each query entry to stdout
    print_best_hits(best_hits_dict)


if __name__ == '__main__':
    main()