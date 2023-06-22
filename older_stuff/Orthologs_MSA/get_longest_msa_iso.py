#!/usr/bin/env python3
#
import argparse
import pyfastx
from operator import itemgetter
from glob import glob
import sys


__author__ = "Ekaterina Osipova, 2020."


def read_fasta_dict(fasta_file):
    ## Reads alignment fasta into a dictionary

    fasta_dict = {}
    for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fasta_dict[name] = seq
    return fasta_dict


def get_msa_length(msa_file, ref):
    ## Reads fasta MSA file; returns its length (by ref)

    fasta_dict = read_fasta_dict(msa_file)
    if ref == '':
        # get any element from the fasta dict
        return len(next(iter(fasta_dict.values())))
    else:
        return len(fasta_dict[ref])


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--msadir', type=str, help='dir with all transcript MSAs sorted by gene')
    parser.add_argument('-r', '--ref', type=str, default='', help='check sequence of this species in each MSA; default: any')
    parser.add_argument('-s', '--suffix', type=str, help='MSA file suffix, e.g: .30birds.prank.hmm.manual.fa')
    args = parser.parse_args()

    ## Make a list of genes in the MSA dir
    msa_dir = args.msadir.rstrip('/')
    gene_dirs_list = [g for g in glob(msa_dir + '/*')]

    ## Make a dictionary: {gene: [list of corresponding MSAs]}
    gene_msa_dict = {g: glob(g + '/*') for g in gene_dirs_list if len(glob(g + '/*')) > 0}

    ## Get MSA size for each transcript for each gene
    ref = args.ref
    gene_msa_size_dict = {}
    for gene in gene_msa_dict:
        gene_msa_size_dict[gene] = [(f, get_msa_length(f, ref)) for f in gene_msa_dict[gene]]

    ## Find the longest MSA for each gene
    suf = args.suffix
    for gene in gene_msa_size_dict:
        msa_size_list = gene_msa_size_dict[gene]
        longest_msa = max(msa_size_list, key=itemgetter(1))[0]
        transcript_name = longest_msa.replace(suf, '').replace(gene + '/', '')
        print('{}\t{}'.format(gene.replace(msa_dir + '/', ''), transcript_name))


if __name__ == "__main__":
        main()
