#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple python script to parse chained homologies from mummer delta file
"""
import os
import sys
import csv
import argparse
import operator
from collections import Counter


def _insert_into_data_struct(name, row, hit_dict):
    """modify dictionary values"""
    if not name in hit_dict:
        hit_dict[name] = [row]
    else:
        hit_dict[name].append(row)
    return hit_dict


def build_dict(matches):
    """build match dictionary"""
    hit_dict = {}
#
    #     for match in matches:
    for match in matches[:7500]:
        rstart, rend, qstart, qend, rlen, qlen, identity, rname, qname = match.split()
        if rname == qname:
            pass
        else:
            row = (qname, rstart, rend, qstart, qend, rlen, qlen, identity)
            hit_dict = _insert_into_data_struct(rname, row, hit_dict)
    return hit_dict


def summary(seen):
    """Generate summary of results"""
    sorted_homologies = [(row[0][0], row[0][1], row[1][0], row[1][1])
                         for row in sorted(seen.items(), key=operator.itemgetter(1))]

    for ref, query, total_bp, chain_length in sorted_homologies:
        print "{x} and {y} share {w} blocks of homology spanning {z} bases".format(x=ref,
                                                                                   y=query,
                                                                                   w=chain_length,
                                                                                   z=total_bp)

    ref_sum_total_bp = sum([i[0] for i in seen.values()])
    print "{x} pairwise(non-reciprocal) shared homologies between contigs spanning {y} bp"\
        .format(x=len(sorted_homologies),
                y=ref_sum_total_bp)

    return sorted_homologies


def write_summary(homologies):
    """write shared_homologies.out
       summary file listing the number of bases / chained homologies for each pair of alignments
    """
    return _write(homologies, ext='out')


def write_coords(homology_list):
    """write shared_homologies.coords
       parsed version of input file containing only homologies of certain chain length and over
    """
    return _write(homology_list, ext='coords')


def _write(items, ext):
    """internal writing function"""
    output = os.path.join(os.getcwd(), 'shared_homologies.{x}'.format(x=ext))

    with open(output, 'wb') as csvfile:
        output = csv.writer(csvfile, delimiter=' ', quotechar=' ')
        for row in items:
            output.writerow(row)


def plot_homologies(seen_dict, coords_file):
    """generate plot commands that can be run afterwards"""
    # writes a bunch of mummerplot commands to a bash script for qsub"""
    path = os.path.dirname(coords_file)
    delta = "{a}{e}".format(a=os.path.basename(coords_file).rstrip('coords'),
                            e='delta')
    delta_path = os.path.join(path, delta)
    cmds = []
    if not os.path.exists(delta_path):
        print 'Delta file not found, cannot generate plots'
    else:

        keys = seen_dict.keys()
        for ref, query in keys:
            output = "{r}_{q}".format(r=ref, q=query)
            cmd = "mummerplot -postscript -r \"{r}\" -q \"{q}\" -p \"{o}\" {d}\n"\
                .format(r=ref,
                        q=query,
                        o=output,
                        d=delta_path)
            cmds.append(cmd)

    output = os.path.join(os.getcwd(), 'plots.sh')
    with open(output, 'wb') as fout:
        fout.writelines(cmds)


def get_homologies(lst, chain_length):
    """get chains of chain_length or higher"""
    ids = [i[0] for i in lst]
    return [i for i in Counter(ids).items() if i[1] >= chain_length]


def get_chain_dict(dictionary):
    """Get dictionary of chain lengths, removing reciprocol alignments"""
    seen = {}
    for reference, homologies in dictionary.iteritems():
        ref_total_bp = 0
        last_query = None

        for alignment in homologies:
            key = (reference, alignment[0])
            if (alignment[0], reference) not in seen.keys():
                if last_query != alignment[0]:
                    chain_len = 1
                    ref_total_bp = int(alignment[5])
                elif last_query == alignment[0]:
                    chain_len += 1
                    ref_total_bp += int(alignment[5])
                seen[key] = (ref_total_bp, chain_len)
            last_query = alignment[0]
    return seen


def parse_homologies(coords, chain_length):
    """Parse chains of homologies"""
    with open(coords, 'r') as fout:
        coords_data = fout.read().splitlines()

    hit_dict = build_dict(coords_data)
    homology_list = []
    for k, value in hit_dict.iteritems():
        homologies = get_homologies(value, chain_length)
        for name, _ in homologies:
            rows = [(k, ) + i for i in [z for z in value if z[0] == name]]
            homology_list.extend(rows)

    sdict = {i[0]: [z[1:] for z in homology_list if z[0] == i[0]]
             for i in homology_list}

    seen = get_chain_dict(sdict)

    return homology_list, seen


def main():
    """run the chain filter, summarize results and generate output"""

    parser = get_parser()
    args = parser.parse_args()
    coords_file = args.coords
    chain_length = args.chain_length

    homology_list, seen = parse_homologies(coords_file, chain_length)
    homologies = summary(seen)
    plot_homologies(seen, coords_file)
    write_summary(homologies)
    write_coords(homology_list)


def get_parser():
    """Parse input"""
    desc = "Tool for parsing chained homologies based on mummer coords file"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('coords', type=str,
                        help="Path to mummer show-coords output")
    parser.add_argument('--chain_length', type=int, default=5,
                        help="Minimum number of chain links")

    return parser

if __name__ == "__main__":

    sys.exit(main())
