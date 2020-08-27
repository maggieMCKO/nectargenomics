#!/usr/bin/env python3
#

import sys
import argparse
from collections import defaultdict
from collections import Counter

__author__ = "Ekaterina Osipova, 2020."


def parse_orthologs(orthologs_file):
    ## Reads orthologs.tsv file into a dictionary {transcript name : [list of orthologous chains]}

    orthologs_dict = {}
    with open(orthologs_file, 'r') as inf:
        for line in inf.readlines():
            # skip the header and empty lines
            if line.startswith('GENE') or not line.strip():
                pass
            else:
                transcript = line.split()[0]
                orthologous_chains = line.split()[1].split(',')
                orthologs_dict[transcript] = orthologous_chains
    return orthologs_dict


def parse_bed_exons(anno_bed_file):
    ## Extracts exon number for each transcript from annotation file

    exon_number_dict = {}
    with open(anno_bed_file, 'r') as inf:
        for anno_line in inf.readlines():
            transcript = anno_line.split()[3]
            exon_number = int(anno_line.split()[9])
            exon_number_dict[transcript] = exon_number
    return exon_number_dict


def parse_meta_data(metadata_file, orthologs_dict, type):
    ## Reads meta_data.tsv file into a dictionary {transcript name : [list of exons]}

    transcript_dict = defaultdict(list)
    with open(metadata_file, 'r') as inf:
        for exon_line in inf.readlines():
            # skip the header and empty lines
            if (exon_line.startswith('gene')) or not (exon_line.strip()):
                pass
            else:
                transcript = exon_line.split()[0]
                exon_num = int(exon_line.split()[1])
                chain_id = exon_line.split()[2]
                predict_class = exon_line.split()[8]
                # check if this is an orthologous chain and the class is what user expects
                if (chain_id in orthologs_dict[transcript]):
                    if (type == 'A') and (predict_class == 'A'):
                        # only A
                        transcript_dict[transcript].append(exon_num)
                    elif (type == 'B') and ((predict_class == 'A') or (predict_class == 'B')):
                        # A and B
                        transcript_dict[transcript].append(exon_num)
                    elif (type == 'M') and ((predict_class == 'A') or (predict_class == 'B') \
                                                    or (predict_class == 'M')):
                        # A, B and M
                        transcript_dict[transcript].append(exon_num)
    return transcript_dict


def classify_by_majority():
    ## Finds the number of fragments the majority exons are covered with;
    ## returns a list of one element (2 elements for ties)

    count_dict = Counter(exon_counts_list)
    max_count_list = [i for i in count_dict if count_dict[i] == max(count_dict.values())]
    return max_count_list


def classify_by_subtract(exon_count_list):
    ## Counts how many times an orthologous reagion can be 'subtracted' from majority of exons
    ## NB: deals with ties conservatively: 2 2 1 1 -> is one-2-one

    count = 0
    while (exon_count_list.count(0) <= round(float(len(exon_count_list))/float(2))) and (exon_count_list != [0]):
        exon_count_list = [i - 1 if i > 0 else 0 for i in exon_count_list]
        count += 1
    return count


def classify_trancripts(transcript_dict, exon_number_dict):
    ## Assigns number of potential orthologs for each transcript
    ## currently uses subtract approach

    for transcript in transcript_dict:
        exon_number = exon_number_dict[transcript]
        #print(exon_number)
        exon_counts_list = [transcript_dict[transcript].count(i) for i in range(exon_number)]
        #print(exon_counts_list)
        transcript_class = classify_by_subtract(exon_counts_list)
        # create output with transcript - class table
        print('{}\t{}'.format(transcript, transcript_class))
    return


def print_exon_counts(transcript_dict, exon_number_dict):
    ## Creates output like: transcript name 1 1 1 2 1 1 2 (exon counts)

    for transcript in transcript_dict:
        exon_number = exon_number_dict[transcript]
        exon_counts_list = [transcript_dict[transcript].count(i) for i in range(exon_number)]
        exon_counts_line = '\t'.join(map(str, exon_counts_list))

        ## to output dictionary structure
        # exon_counts_dict = Counter(transcript_dict[transcript])
        # exon_counts_line = '\t'.join(map(str, [exon_counts_dict[i] for i in exon_counts_dict]))
        print('{}\t{}'.format(transcript, exon_counts_line))
    return


def main():
    ## Parse arguments
    parse = argparse.ArgumentParser()
    parse.add_argument('-o', '--orthologs', type=str, help='orthologs.tsv file produced by TOGA')
    parse.add_argument('-m', '--metadata', type=str, help='meta_data.tsv file produced by TOGA')
    parse.add_argument('-b', '--bed12', type=str, help='annotation file in bed12 format')
    parse.add_argument('-t', '--type', type=str, default='A', choices={'A', 'B', 'M'}, help='class, A for only A, B - for A and B; M - for A,B and M')
    args = parse.parse_args()

    ## Read orthologs file into a dictionary
    orthologs_dict = parse_orthologs(args.orthologs)

    ## Read annotation bed file into a dictionary: {transcript : number of exons}
    exon_number_dict = parse_bed_exons(args.bed12)

    ## Read meta data file into a dictionary
    transcript_dict = parse_meta_data(args.metadata, orthologs_dict, type=args.type)

    ## Output counts of each exon of each transcript
    #print_exon_counts(transcript_dict, exon_number_dict)

    ## Classify transcripts using counts of its exons; output it as a table
    classify_trancripts(transcript_dict, exon_number_dict)


if __name__ == '__main__':
    main()







