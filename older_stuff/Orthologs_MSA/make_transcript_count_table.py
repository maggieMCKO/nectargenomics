#!/usr/bin/env python3
#
import argparse
from collections import defaultdict
import re
from os import listdir
from os.path import isfile, join


__author__ = "Ekaterina Osipova, 2020."



def read_species_from_file(species_file):
    ## Read species file into a list

    species_list = []
    with open(species_file, 'r') as inf:
        for line in inf.readlines():
            species_list.append(line.rstrip())
    return list(set(species_list))


def read_counts_files_into_dict(counts_dir):
    ## Reads all counts files in the counts dir into a trasnscript dictionary

    transcript_dict = defaultdict(list)
    for count_file in [join(counts_dir, f) for f in listdir(counts_dir) if isfile(join(counts_dir, f))]:

        # get db name
        db_pattern = 'H?L?[a-z]{3}[A-Z]{1}[a-z]{2}[0-9]{1}'
        m = re.search(db_pattern, count_file)
        db = m.group(0)

        # read file into dictionary
        with open(count_file, 'r') as inf:
            for line in inf.readlines():
                transcript, count = line.split()[0], line.split()[1]
                transcript_dict[transcript].append((db, count))
    return transcript_dict


def convert_dict_table(transcript_dict, species_list):
    ## Makes a table with rows with transcript ids and columns with counts of a transcript for each species

    print('\t' + '\t'.join(species_list))
    for transcript in transcript_dict:
        transcript_line = '{}\t'.format(transcript)
        for db in species_list:
            if db in [db_count[0] for db_count in transcript_dict[transcript]]:
                db_index = [db_count[0] for db_count in transcript_dict[transcript]].index(db)
                count = transcript_dict[transcript][db_index][1]
            else:
                count = '0'
            transcript_line += '{}\t'.format(count)
        print(transcript_line)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--specieslist', type=str, help='file with a list of species')
    parser.add_argument('-c', '--countsdir', type=str, help='dir with exon.counts.$db.class.tsv files')
    args = parser.parse_args()

    ## Read species file into a list
    species_list = read_species_from_file(args.specieslist)

    ## For each species read the exon.counts file
    transcript_dict = read_counts_files_into_dict(args.countsdir)

    ## Prepare a table with counts for each transcript for each species
    convert_dict_table(transcript_dict, species_list)


if __name__ == "__main__":
    main()