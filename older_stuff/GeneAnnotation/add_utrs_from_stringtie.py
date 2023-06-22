#!/usr/bin/env python
#

"""
This script takes annotation in bed12 format without UTRs and adds UTRs from a transcriptome assembly (e.g. stringtie)
"""

import argparse
from collections import defaultdict
from collections import Counter
import sys

__author__ = "Ekaterina Osipova, 2020."


def get_all_intron_coords(transcript):
    ## Takes bed12 transcript line and returns a list with intron intervals: [(chrom, x1, x2), (chrom, x3, x4), ..]

    chrom = transcript.split()[0]
    start = int(transcript.split()[1])
    exon_number = int(transcript.split()[9])
    intron_coord_list = []
    block_starts = [start + i for i in map(int, transcript.split()[11].rstrip(',').split(','))]
    block_sizes = map(int, transcript.split()[10].rstrip(',').split(','))
    for i in range(exon_number - 1):
        intron_coord_list.append((chrom, block_starts[i] + block_sizes[i], block_starts[i + 1]))
    return intron_coord_list


def read_rnaseq_bed(rnaseq_file):
    ## Reads RNAseq annotation file (no CDS regions) and adds all introns to rnaseq_coord_dict:
    ## {(chrom, x1, x2): [transcript1, transcript2, ..], }

    rnaseq_coord_dict = defaultdict(list)
    with open(rnaseq_file, 'r') as inf:
        for bed_line in inf.readlines():
            intron_coord_list = get_all_intron_coords(bed_line)
            for intron in intron_coord_list:
                rnaseq_coord_dict[intron].append(bed_line)
    return rnaseq_coord_dict


def prepare_5prime_blocks(transcript, cut_coord):
    ## Takes bed12 transcript line and one coord where to cut it;
    ## returns list of starts and sizes of non-coding blocks; works with 5'-UTR

    noncoding_blocks = []
    start = int(transcript.split()[1])
    abs_block_starts_list = [start + int(i) for i in transcript.split()[11].rstrip(',').split(',')]
    block_sizes_list = [int(i) for i in transcript.split()[10].rstrip(',').split(',')]

    i = 0
    while abs_block_starts_list[i] + block_sizes_list[i] < cut_coord:
        noncoding_blocks.append((abs_block_starts_list[i], block_sizes_list[i]))
        i += 1
    # add the cut block
    noncoding_blocks.append((abs_block_starts_list[i], cut_coord - abs_block_starts_list[i]))
    return noncoding_blocks


def prepare_3prime_blocks(transcript, cut_coord):
    ## Takes bed12 transcript line and one coord where to cut it;
    ## returns list of starts and sizes of non-coding blocks; works with 3'-UTR

    noncoding_blocks = []
    start = int(transcript.split()[1])
    abs_block_starts_list = [start + int(i) for i in transcript.split()[11].rstrip(',').split(',')]
    block_sizes_list = [int(i) for i in transcript.split()[10].rstrip(',').split(',')]

    abs_block_starts_list.reverse()
    block_sizes_list.reverse()
    i = 0
    while abs_block_starts_list[i] > cut_coord:
        noncoding_blocks.append((abs_block_starts_list[i], block_sizes_list[i]))
        i += 1
    # add the cut block
    noncoding_blocks.append((cut_coord, block_sizes_list[i] - cut_coord + abs_block_starts_list[i]))
    # do you need to reverse it back??
    noncoding_blocks.reverse()
    return noncoding_blocks


def read_anno_bed(anno_file, rnaseq_coord_dict):
    ## Reads annotation file without UTRs and finds if first/last intron matches perfectly anything in rnaseq_coord_dict

    transcript_dict = defaultdict(list)
    with open(anno_file, 'r') as inf:
        for bed_line in inf.readlines():
            name = bed_line.split()[3]
            start = int(bed_line.split()[1])
            end = int(bed_line.split()[2])
            exon_number = int(bed_line.split()[9])
            transcript_info = bed_line

            ## Initiate and update UTR lists in the transcript_dict (3'UTRs and 5'UTRs separately)
            utrs5_list = []
            utrs3_list = []

            if exon_number > 1:
                intron_coord_list = get_all_intron_coords(bed_line)
                first_intron_coords = intron_coord_list[0]
                last_intron_coords = intron_coord_list[-1]

                if (first_intron_coords in rnaseq_coord_dict):
                    utrs5_list = [i for i in rnaseq_coord_dict[first_intron_coords] if int(i.split()[1]) < start]
                if (last_intron_coords in rnaseq_coord_dict):
                    utrs3_list = [i for i in rnaseq_coord_dict[last_intron_coords] if int(i.split()[2]) > end]

            transcript_dict[name].append((transcript_info, utrs5_list, utrs3_list))
    return transcript_dict


def get_max_most_common(list, utr_type):
    ## Finds most common element in the list; if there're multiple, returns max of those

    count_dict = Counter(list)
    maxcount = count_dict.most_common(1)[0][1]
    if utr_type == 5:
        max_most_common = min([i for i in count_dict if count_dict[i] == maxcount])
    else:
        max_most_common = max([i for i in count_dict if count_dict[i] == maxcount])
    return max_most_common



def add_utr_blocks(transcript_info, utr5_blocks, utr3_blocks):
    ## Given transcript annotation line updates it with utrs

    start = int(transcript_info.split()[1])
    chrom = transcript_info.split()[0]
    exon_number = int(transcript_info.split()[9])
    name = transcript_info.split()[3]
    trans_starts = [i + start for i in map(int, transcript_info.split()[11].rstrip(',').split(','))]
    trans_sizes = map(int, transcript_info.split()[10].rstrip(',').split(','))

    utr5_starts = [i[0] for i in utr5_blocks]
    utr5_sizes = [i[1] for i in utr5_blocks]
    utr3_starts = [i[0] for i in utr3_blocks]
    utr3_sizes = [i[1] for i in utr3_blocks]


    start_update = utr5_starts[0]
    end_update = utr3_starts[-1] + utr3_sizes[-1]

    ## update starts and block sizes
    new_trans_starts = [i - utr5_starts[0] for i in utr5_starts + trans_starts[1:] + utr3_starts[1:]]
    block_starts_update = ','.join(map(str, new_trans_starts)) + ','

    # make updates of block sizes if it's a single exon gene
    if exon_number == 1:
        new_trans_sizes = utr5_sizes[:-1] + [utr5_sizes[-1] + trans_sizes[0] + utr3_sizes[0]] + utr3_sizes[1:]
    # make updates of block sizes if number of exons > 1
    else:
        new_trans_sizes = utr5_sizes[:-1] + [utr5_sizes[-1] + trans_sizes[0]] + trans_sizes[1:-1] + \
                            [trans_sizes[-1] + utr3_sizes[0]] + utr3_sizes[1:]
    block_sizes_update = ','.join(map(str, new_trans_sizes)) + ','

    ## combine all info in a new annotation line
    exon_number_update = len(new_trans_starts)
    bed_line_list_update = [chrom, str(start_update), str(end_update), name] + transcript_info.split()[4:9] + \
                           [str(exon_number_update), block_sizes_update, block_starts_update]
    bed_line_update = '\t'.join(bed_line_list_update)

    return bed_line_update


def update_annotation(transcript_dict):
    ## Adds 5'- and 3'-UTRs for each transcript in given transcript_dict
    ## Runs add_utrs() function that work with an individual transcript

    for name in transcript_dict:
        for transcript in transcript_dict[name]:
            transcript_info = transcript[0]
            cut_utr5 = int(transcript_info.split()[1])
            cut_utr3 = int(transcript_info.split()[2])
            utr5_transcript_list = transcript[1]
            utr3_transcript_list = transcript[2]

            # check if all non-coding starts > start !!! not yet done
            # print('utr5 list: ', utr5_transcript_list)
            # print('utr3 list: ', utr3_transcript_list)

            # if there are updates for this transcript, get the most common coordinate
            if utr5_transcript_list != []:
                utr5_start_list = [int(i.split()[1]) for i in utr5_transcript_list]
                utr5 = get_max_most_common(utr5_start_list, utr_type=5)
                maxcount_transcript_index = [int(i.split()[1]) for i in utr5_transcript_list].index(utr5)
                utr5_blocks = prepare_5prime_blocks(utr5_transcript_list[maxcount_transcript_index], cut_utr5)
            else:
                utr5_blocks = [(int(transcript_info.split()[1]), 0)]

            if utr3_transcript_list != []:
                utr3_end_list = [int(i.split()[2]) for i in utr3_transcript_list]
                utr3 = get_max_most_common(utr3_end_list, utr_type=3)
                maxcount_transcript_index = [int(i.split()[2]) for i in utr3_transcript_list].index(utr3)
                utr3_blocks = prepare_3prime_blocks(utr3_transcript_list[maxcount_transcript_index], cut_utr3)
            else:
                utr3_blocks = [(int(transcript_info.split()[2]), 0)]

            bed_line_update  = add_utr_blocks(transcript_info, utr5_blocks, utr3_blocks)
            print(bed_line_update)
    return


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--rnaseq', type=str, help='transcripts assembled from RNAseq, e.g. stringtie in bed12 format')
    parser.add_argument('-a', '--anno', type=str, help='bed12 annotation file to add UTRs to')
    args = parser.parse_args()

    ## bed12 format:
    ## chrom[0] start[1] end[2] name[3] score[4] strand[5] cds_start[6] cds_end[7] rgb[8] count[9]\
    ##  block_sizes[10] block_starts[11]

    ## Read stringtie assembly into a dictionary
    rnaseq_coord_dict = read_rnaseq_bed(args.rnaseq)

    ## Read annotation file checking if first/last blocks overlap blocks in rnaseq
    transcript_dict = read_anno_bed(args.anno, rnaseq_coord_dict)

    ## Add UTRs to the transcripts where possible
    update_annotation(transcript_dict)


if __name__ == "__main__":
    main()