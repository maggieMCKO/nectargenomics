#!/usr/bin/env python
#
# This script takes two files with set of elements in format element \n \
# and outputs what in not in the first set but in the second; what is not in the second but in the first; \
# length of these differences between two sets.
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f1', '--fileone', type=str, help='file with first set of elements')
parser.add_argument('-f2', '--filetwo', type=str, help='file with second set of elements')

args = parser.parse_args()

with open(args.fileone, 'r') as fone:
    lines = fone.readlines()
    firstList = map(str.strip, lines)

with open(args.filetwo, 'r') as ftwo:
    lines = ftwo.readlines()
    secondList = map(str.strip, lines)

seen = set()
seenAdd = seen.add
uniqFirst = [x for x in secondList if not (x in seen or seenAdd(x))]
seen = set()
seenAdd = seen.add
uniqSecond = [x for x in firstList if not (x in seen or seenAdd(x))]

gainList = list(set(uniqSecond) - set(uniqFirst))
lossList = list(set(uniqFirst) - set(uniqSecond))

print '# of elements NOT in f2: ', len(list(set(uniqSecond) - set(uniqFirst)))
print '# of elements NOT in f1: ', len(list(set(uniqFirst) - set(uniqSecond)))
print '\n'
print 'Elements NOT in f2: \n', gainList
print '\n'
print 'Elements NOT in f1: \n', lossList
