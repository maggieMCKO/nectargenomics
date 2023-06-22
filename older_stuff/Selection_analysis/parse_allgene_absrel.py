#!/usr/bin/env python3
#
'''
    This script parses all.gene absrel table in format branch \t pval \t transcript \t gene
    It can calculate difference between p-values for transcripts fo each gene;
    It can output transcript with the lowest p-value for each gene
'''

import argparse
from collections import defaultdict
import sys
import itertools


__author__ = "Ekaterina Osipova, 2020."


def parse_pvalue_table(file):
    ## Reads input all gene table

    gene_branch_dict = defaultdict(list)
    gene_dict = defaultdict(list)
    with open(file, 'r') as inf:
        for line in inf.readlines():
            if not (line.startswith('branch')):
                branch = line.split()[0]
                pval = float(line.split()[1])
                trans = line.split()[2]
                gene = line.split()[3]
                gene_branch_dict[(gene, branch)].append((trans, pval))
                gene_dict[gene].append((branch, trans, pval))
    return gene_branch_dict, gene_dict


def calc_diff(f_list):
    ## Calculates difference between each pair of elements if the float list

    diff_list = []
    for pair in itertools.combinations(f_list, 2):
        diff_list.append(abs(pair[0] - pair[1]))
    return diff_list


def get_genes_diff_pvalues(gene_branch_dict, delta, pcutoff):
    ## Outputs genes that are having different p-values for different transcripts

    for gene_branch in gene_branch_dict:
        transc_pval_list = gene_branch_dict[gene_branch]
        pval_list = [i[1] for i in transc_pval_list]
        if len(pval_list) > 1:
            pval_diff_list = calc_diff(pval_list)
            bool_pval_diff_list = [True if i >= delta else False for i in pval_diff_list]
            bool_pval_cutoff_list = [True if i <= pcutoff else False for i in pval_list]
            if any(bool_pval_cutoff_list) and any(bool_pval_diff_list):
                transc_pval_str = ';'.join([i[0] + ':' + str(i[1]) for i in transc_pval_list])
                print('{}\t{}\t{}'.format(gene_branch[1], transc_pval_str, gene_branch[0]))


def get_most_extreme_transc(gene_dict, mode):
    ## Finds the lowest p-value transcript for each gene; returns in form of dictionary

    extreme_gene_dict = {}
    for gene in gene_dict:
        pval_list = [i[2] for i in gene_dict[gene]]
        if len(pval_list) > 1:
            if mode == 'lowest':
                extreme_pval = min(pval_list)
            elif mode == 'highest':
                extreme_pval = max(pval_list)
            else:
                print('Provide appropriate mode!: lowest/highest/diff')
                sys.exit(1)
        else:
            extreme_pval = pval_list[0]
        extreme_transc = [i[1] for i in gene_dict[gene] if i[2] == extreme_pval][0]
        extreme_gene_dict[gene] = extreme_transc
    return extreme_gene_dict


def print_most_extreme(extreme_gene_dict, gene_branch_dict, mode):
    ## If there are inconsistencies in the transcripts choice, get the lowest of the LOWEST

    for gene_branch in gene_branch_dict:
        gene = gene_branch[0]
        branch = gene_branch[1]
        extreme_transc =  extreme_gene_dict[gene]
        if extreme_transc in [i[0] for i in  gene_branch_dict[gene_branch]]:
            extreme_pval = [i[1] for i in gene_branch_dict[gene_branch] if i[0] == extreme_transc][0]
            print('{}\t{}\t{}\t{}'.format(branch, extreme_pval, extreme_transc, gene))


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_alltab', type=str, help='all gene table: branch\tpval\ttrans\tgene')
    parser.add_argument('-m', '--mode', type=str, default='lowest', help='mode you want to run the scripts in; \
                        lowest/highest: to get lowest/highest p-value; diff: to get p-value differences')
    ## Arguments for p-value difference mode
    parser.add_argument('-d', '--delta', type=float, default=.01, help='min required difference in p-values; default: .01')
    parser.add_argument('-p', '--pcutoff', type=float, default=.05, help='at least one of the transcripts for a gene \
                        should have lower p-value')
    args = parser.parse_args()

    ## Parse input table into a dictionary
    gene_branch_dict, gene_dict = parse_pvalue_table(args.input_alltab)

    ## Check the mode to run the script in and run it
    mode = args.mode
    if mode == 'lowest' or mode == 'highest':
        ## Print a table in the input format with only highest/lowest p-value transcript for each gene
        extreme_gene_dict = get_most_extreme_transc(gene_dict, mode)
        print_most_extreme(extreme_gene_dict, gene_branch_dict, mode)
    else:
        ## Output genes that are having different p-values for different transcripts
        get_genes_diff_pvalues(gene_branch_dict, args.delta, args.pcutoff)


if __name__ == "__main__":
        main()