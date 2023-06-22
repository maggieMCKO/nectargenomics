#!/usr/bin/env python
#
# This script takes files with set of elements in format element \n \
# and outputs the union of elements; \
import argparse


__author__ = "Ekaterina Osipova, 2020."


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filelist', nargs='*', type=str, help='list of files to find union')
args = parser.parse_args()


file_list = args.filelist
allElements_allFiles = []

for file in file_list:
    with open(file, 'r') as inf:
        for line in inf.readlines():
            allElements_allFiles.append(line.strip())

union_elements = set(allElements_allFiles)

# output union of elements
for el in union_elements:
    print(el)
