#!/usr/bin/env python3
#
import argparse
import json


__author__ = "Ekaterina Osipova, 2020."


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--absreljson', type=str, help='file with json output of aBSREL')
    parser.add_argument('-a', '--attribute', type=str, default='Corrected P-value', help='attribute of a branch you want to extract; default: Corrected P-value')
    parser.add_argument('-b', '--branches', type=str, default='', help='output only these branches (comma-sep); default: all')
    args = parser.parse_args()

    ## Read aBSREL output into json object
    with open(args.absreljson, 'r') as inf:
        data = json.load(inf)

    ## Get requested attribute
    attribute = args.attribute
    branch_attr = data['branch attributes']['0']
    out_dict_branch_attr = {}
    for branch in branch_attr:
        attr_value = branch_attr[branch][attribute]
        out_dict_branch_attr[branch] = attr_value

    ## Output attribute for requested branches
    if args.branches == '':
        for branch in out_dict_branch_attr:
            attr_value = out_dict_branch_attr[branch]
            print('{}\t{}'.format(branch, attr_value))
    else:
        for branch in args.branches.split(','):
            attr_value = out_dict_branch_attr[branch]
            print('{}\t{}'.format(branch, attr_value))


if __name__ == "__main__":
    main()
