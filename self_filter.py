#!/usr/bin/env python
"""When identifying contig -> contig relationships in a single assembly
   fasta, this simple script removes self hits from the delta file
   allowing for more effective `delta-filter`ing

"""

import sys
import argparse


def filter_exact_delta(delta_data, output):
    """remove self hits"""
    new_delta = []
    for i, row in enumerate(delta_data):
        if i == 0 or i == 1:
            new_delta.append(row)
        else:
            if row.startswith(">"):
                self_hit = False
                ref, query, _, _ = row.lstrip('>').split()
                if ref == query:
                    self_hit = True
                else:
                    new_delta.append(row)

            else:

                if self_hit is True:
                    pass
                else:
                    new_delta.append(row)

    write_output(new_delta, output)

    return new_delta


def write_output(new_delta, output):
    """write noself delta"""
    output = "{x}_noself.delta".format(x=output.rstrip('.delta'))
    with open(output, 'w') as fout:
        for item in new_delta:
            fout.write("{i}\n".format(i=item))


def main():
    """run self hit filter"""
    parser = get_parser()
    args = parser.parse_args()
    delta = args.delta
    output = "{x}_noself.delta".format(x=delta.rstrip('.delta'))

    with open(delta, 'r') as fin:
        rows = fin.read().splitlines()

    filter_exact_delta(rows, output)

    return


def get_parser():
    """Parse input"""
    desc = "Simple python script to remove self hits from a delta file " \
           "when aligning a set of contigs back to itself"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('delta',
                        type=str,
                        help="Path to mummer delta-file containing self alignments")

    return parser


if __name__ == "__main__":
    sys.exit(main())
