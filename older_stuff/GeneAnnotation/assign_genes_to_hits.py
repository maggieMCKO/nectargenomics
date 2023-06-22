#!/usr/bin/env python3
#

"""
This script substitutes uniprot IDs in a blast hits file with corresponding gene names (where possible!).
example output: ENSGALT00000098364.1,BCL9L
"""

import argparse

__author__ = "Ekaterina Osipova, 2020."


def read_uniprot_db(uniprot):
    ## Reads uniprot db into a dictionary

    uniprot_db_dict = {}
    with open(uniprot, 'r') as inf:
        for line in inf.readlines():
            if line.startswith('>'):
                uniprot_ids = [i for i in line.rstrip().split() if i.startswith('>sp')]
                uniprot_id = uniprot_ids[0].split('|')[1]
                gene_names = [i for i in line.rstrip().split() if i.startswith('GN=')]
                if len(gene_names) > 0:
                    gene_name = gene_names[0].lstrip('GN=')
                    uniprot_db_dict[uniprot_id] = gene_name
                else:
                    uniprot_db_dict[uniprot_id] = 'NO_UNIREF_GENE'
    return uniprot_db_dict


def assign_genes_to_hits(hits_file, uniprot_db_dict, suffix):
    ## Reads blast hits and assigns gene names of corresponding uniprot IDs to transcripts

    with open(hits_file, 'r') as inf:
        for line in inf.readlines():
            transcript = line.split()[0]
            uniprot_id = line.split()[1].split('|')[1]
            if uniprot_db_dict[uniprot_id] != 'NO_UNIREF_GENE':
                print('{},{}-{}'.format(transcript, uniprot_db_dict[uniprot_id], suffix))
    return

def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--uniprotdb', type=str, help='file with uniprot/swissprot database')
    parser.add_argument('-b', '--blasthits', type=str, help='file with blast hits in outfmt 6')
    parser.add_argument('-s', '--suffix', type=str, default='', help='suffix for gene names you can add, e.g: FGL1-like')
    args = parser.parse_args()

    ## Read uniprot/swissprot database into a dictionary {ID: gene_name}
    uniprot_db_dict = read_uniprot_db(args.uniprotdb)

    ## Make csv table of transcripts and corresponding gene names; stdout
    assign_genes_to_hits(args.blasthits, uniprot_db_dict, args.suffix)


if __name__ == '__main__':
    main()