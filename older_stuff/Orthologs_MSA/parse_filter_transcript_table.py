#!/usr/bin/env python3
#
import argparse
import re



__author__ = "Ekaterina Osipova, 2020."


def parse_transcript_table(table_counts, required_species, fraction, maxcount):

    with open(table_counts, 'r') as inf:
        for line in inf.readlines():
            if re.match(r'\s', line):
                # read header line
                species_line = line.strip().split()
                species_number = len(species_line)
                if required_species != ['NA']:
                    required_species_indexes = [species_line.index(db) for db in required_species]
            else:
                # read transcript line
                transcript, count_list = line.split()[0], [int(i) for i in line.split()[1:]]

                # check if the line doesn't have counts >maxcount and at least $fraction of species has a count >=1
                if (set(count_list).issubset(set(range(maxcount + 1)))) and (count_list.count(0) < (fraction * species_number)):
                    # check if all required species have count >=1
                    if (required_species != ['NA']):
                        if (set([count_list[ind] for ind in required_species_indexes]).issubset(set(range(1, maxcount + 1)))):
                            print('{}'.format(transcript))
                    else:
                        print('{}'.format(transcript))


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tablecount', type=str, help='table of transcripts counts')
    parser.add_argument('-s', '--species', type=str, default='NA', help='list of species required to have count >=1 (comma-sep)')
    parser.add_argument('-f', '--fraction', type=float, default=0.8, help='fraction of species required to have count >=1')
    parser.add_argument('-m', '--maxcount', type=int, default=1, help='max count that any species can have; default = 1')
    args = parser.parse_args()

    ## Parse transcript table; output transcripts that satisfy conditions
    required_species_list = args.species.split(',')
    parse_transcript_table(args.tablecount, required_species_list, args.fraction, args.maxcount)


if __name__ == "__main__":
    main()