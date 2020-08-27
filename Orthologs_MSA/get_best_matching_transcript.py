#!/usr/bin/env python3
#
import argparse
import numpy as np


__author__ = "Ekaterina Osipova, 2020."


def get_bases(starts, sizes):
    ## Converts provided intervals into lists of bases

    list_2d = [list(range(starts[i], starts[i] + sizes[i])) for i in range(len(starts))]
    return list(np.concatenate(list_2d))


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--transcriptlist', type=str, help='list of transcripts: a,b from bedtools intersect')
    args = parser.parse_args()

    ## Read each bed12+bed12 line into transcript dictionary with coverage of A and B
    cov_dict = {}

    ## bed12 format:
    ## chrom[0] start[1] end[2] name[3] score[4] strand[5] cds_start[6] cds_end[7] rgb[8] count[9]\
    ##  block_sizes[10] block_starts[11]
    with open(args.transcriptlist, 'r') as inf:
        for line in inf.readlines():
            a_bed = line.split()[:12]
            b_bed = line.split()[12:25]
            trans_name = a_bed[3]

            a_starts = [int(i) for i in a_bed[11].rstrip(',').split(',')]
            b_starts = [int(i) for i in b_bed[11].rstrip(',').split(',')]
            a_block_sizes = [int(i) for i in a_bed[10].rstrip(',').split(',')]
            b_block_sizes = [int(i) for i in b_bed[10].rstrip(',').split(',')]
            a_bases = get_bases(a_starts, a_block_sizes)
            b_bases = get_bases(b_starts, b_block_sizes)

            bases_overlap = len(list(set(a_bases) & set(b_bases)))
            a_delta = 1.0 - float(bases_overlap) / float(len(a_bases))
            b_delta = 1.0 - float(bases_overlap) / float(len(b_bases))
            fine = a_delta + b_delta

            #print('A: {}\tB: {}\tfine: {}'.format(a_bed[3], b_bed[3], fine))
            # { rna-XM_030279135.1.4 : fine }
            cov_dict[fine] = trans_name

    ## Find and output the transcript with the best overlap
    best_overlap = cov_dict[min([k for k in cov_dict])]
    print('{}'.format(best_overlap))


if __name__ == "__main__":
    main()