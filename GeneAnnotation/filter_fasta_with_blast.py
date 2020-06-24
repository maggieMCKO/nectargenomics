#!/usr/bin/env python
#

"""
This script takes fasta file with predicted gene models, blast hits file and \
a file with list of species you want to include;
it processes blast hilts and filters out all proteins in fasta file that don't have high-scoring hits\
in blast file or coming from other species to species from the provided list
"""


import argparse
from Bio import SeqIO
from collections import defaultdict
import sys

__author__ = "Ekaterina Osipova, 2020."


def make_allowed_species_list(species_file):
    ## Reads species file and makes a list of allowed species codes: [CHICK, URILO, ..]

    species_codes = []
    with open(species_file, 'r') as inf:
        for line in inf.readlines():
            species_codes.append(line.split('\t')[0])
    return species_codes


def read_blast_hits(blast_file, species_codes, idmin, qcov, rcov):
    ## Reads blast outfmt6 hits into a dictionary if hit passes quality check
    ## qseqid, rseqid, pid, alilen, mism, gapop, qst, qend, rst, rend, eval, bitscore

    blast_dict = defaultdict(list)
    with open(blast_file, 'r') as blastOut:
        for line in blastOut.readlines():
            qseqid = line.split('\t')[0]
            pid = float(line.split('\t')[2])
            ali_length = float(line.split('\t')[3])
            qEnd = float(line.split('\t')[7])
            rEnd = float(line.split('\t')[9])

            ref_code = (line.split('\t')[1]).split('_')[-1]
            if ref_code in species_codes:
                if (pid >= idmin) and (ali_length >= qcov * qEnd) and (ali_length >= rcov * rEnd):
                    entry_info = '\t'.join(line.split('\t')[1:])
                    blast_dict[qseqid].append(entry_info)
    return blast_dict


def filter_fasta(fasta_file, blast_dict, lenmin):
    ## Goes through fasta file and checks each entry for presence in the blast_dictionary (later: for quality!)

    sequences = SeqIO.parse(open(fasta_file), 'fasta')
    for entry in sequences:
        header, seq = entry.id, str(entry.seq)
        if (header in blast_dict) or (len(seq) >= lenmin):
            print('>{}'.format(header))
            print(seq)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--blast', type=str, help='blast output in outfmt 6 format')
    parser.add_argument('-f', '--fasta', type=str, help='fasta file with predicted gene models')
    parser.add_argument('-s', '--species', type=str, help='species file, latin names')
    parser.add_argument('-l', '--lenmin', type=int, help='min required length of a predicted protein if not found in db')
    parser.add_argument('-id', '--idmin', type=float, default=50, help='min required % of identity between query and target')
    parser.add_argument('-qcov', '--qcov', type=float, default=0.75, help='min required query coverage in a hit')
    parser.add_argument('-rcov', '--rcov', type=float, default=0.5, help='min required target coverage in a hit')
    args = parser.parse_args()

    ## Read species file and make a list of allowed species codes: [CHICK, URILO, ..]
    species_codes = make_allowed_species_list(args.species)

    ## Read blast outfmt6 file into a dictionary
    blast_dict = read_blast_hits(args.blast, species_codes, args.idmin, args.qcov, args.rcov)

    ## Check each fasta if it has a blast hit (and quality); output to stdout
    filter_fasta(args.fasta, blast_dict, args.lenmin)


if __name__ == '__main__':
    main()



