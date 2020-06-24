#!/usr/bin/env python3

# this script replaces name of scaffold/chroms in the fist column of the annotation file (format does not matter)
# with corresponding scaffold/chroms names from the renaming_dictionary.csv file
# writes to stdout

import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-a', '--annofile', type=str, help='annotation file; 1st columns - scaffolds/chroms names')
parser.add_argument('-d', '--renamedict', type=str, help='table of name correspondence, usually renaming_dictionary.csv')
args = parser.parse_args()

# read correspondence table into dictionary
renameTable = {}
with open(args.renamedict, 'r') as inf:
    for line in inf.readlines():
        renameTable[line.split(',')[0]] = line.rstrip().split(',')[1]


# read annotation file, replace 1st column
with open(args.annofile, 'r') as inf:
    for line in inf.readlines():
        if line.split('\t')[0] in renameTable:
            newScaffName = renameTable[line.split('\t')[0]]
        else:
            newScaffName = line.split('\t')[0]
        newLine = line.rstrip().split('\t')[1:]
        newLine.insert(0, newScaffName)
        print(*newLine, sep='\t')


