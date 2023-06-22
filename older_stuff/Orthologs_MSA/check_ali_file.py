#!/usr/bin/env python3
#
import argparse
import pyfastx
import re



__author__ = "Ekaterina Osipova, 2020."


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, help='MSA fasta file to check')
    parser.add_argument('-n', '--number', type=int, help='min number of species required to be present in MSA')
    parser.add_argument('-r', '--required', type=str, default='', help='species that need to be in the MSA, comma-sep')
    parser.add_argument('-t', '--text', type=str, help='text to print if fasta passed quality control')
    args = parser.parse_args()

    ## Reads alignment fasta; checks number of species; checks required species
    fasta_dict = {}
    for name, seq in pyfastx.Fasta(args.fasta, build_index=False):
        fasta_dict[name] = seq

    codon_control = set([len(fasta_dict[name]) % 3 for name in fasta_dict])
    required_species = args.required.split(',')

    if (len(fasta_dict) >= args.number) and (codon_control == {0}) and (all(i in fasta_dict for i in required_species)):
        print(args.text)


if __name__ == "__main__":
        main()