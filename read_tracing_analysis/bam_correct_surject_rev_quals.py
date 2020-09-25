#!/usr/bin/env python2.7

"""
Charles Markello
2/12/2019

Description:
    

"""

from __future__ import print_function
import argparse
import sys
import pysam
from collections import defaultdict
import pdb

def parse_args():
    """
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """

    parser = argparse.ArgumentParser(prog='read-tracer', description= __doc__)
    parser.add_argument('-i', '--in_bam', type=str,
        help="Input BAM file to have reverse-strand-aligned reads qual strings reversed")
    parser.add_argument('-o','--out_fixed_bam', type=str, 
        help="Output BAM file containing the reversed qual strings for reverse-strand-aligned reads.")
    
    options = parser.parse_args()

    return options

def main(args):
    
    options = parse_args()

    in_bamfile = pysam.AlignmentFile(options.in_bam, "rb")
    out_bamfile = pysam.AlignmentFile(options.out_fixed_bam, "wb", template=in_bamfile)
    
    for read_record in in_bamfile:
        if read_record.is_reverse:
            new_query_qualities = read_record.query_qualities[::-1]
            read_record.query_qualities = new_query_qualities
            out_bamfile.write(read_record)
        else:
            out_bamfile.write(read_record)
    
    out_bamfile.close()
    in_bamfile.close()
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))


