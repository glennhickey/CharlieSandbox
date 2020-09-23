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
    parser.add_argument('-a', '--in_bamFile1', type=str,
        help="Input test BAM file sorted by read name via 'samtools -n'")
    parser.add_argument('-b', '--in_bamFile2', type=str,
        help="Input baseline BAM file sorted by read name via 'samtools -n'")
    parser.add_argument('-o','--out_filtered_moved_pos_BAM', type=str, 
        help="Output BAM file that will contain reads from the test BAM that satisfy position movement and map quality difference criterias.")
    parser.add_argument('-i','--out_filtered_initial_pos_BAM', type=str, 
        help="Output BAM file that will contain reads from the baseline BAM that satisfy position movement and map quality difference criterias.")
    parser.add_argument('-p','--out_read_pos_diff', type=argparse.FileType('wb'), default=sys.stdout,
        help= "Output tab-delimited list of records in the form of: \
             READNAME'\t'CHROMOSOME POSITION OF READ IN BAMFILE 1'\t'BAMFILE 1 READ QUAL SCORE'\t'CHROMOSOME POSITION OF READ IN BAMFILE 2'\t'BAMFILE 2 READ QUAL SCORE.")
    parser.add_argument('-d','--pos_diff_thresh', type=int, default=150,
        help= "Threshold of read position difference of a read between two BAM files.")
    parser.add_argument('-m','--mapq_diff_thresh', type=int, default=0,
        help= "Threshold of map quality difference of a read between two BAM files.")
    
    options = parser.parse_args()

    return options

def main(args):
    
    options = parse_args()

    # Cycle through read records in each BAM file and extract reads that share the same name and which pair that read is
    #   and filter for reads that have position that either is in a different reference contig or a different position
    #   a different position that's >= 'diff_thresh' within the same reference contig
    print('loading bam files...')
    bamfile_1 = pysam.AlignmentFile(options.in_bamFile1, "r")
    bamfile_2 = pysam.AlignmentFile(options.in_bamFile2, "r")
    bamfile_diff_filtered = pysam.AlignmentFile(options.out_filtered_moved_pos_BAM, "w", template=bamfile_1)
    bamfile_inital_pos_filtered = pysam.AlignmentFile(options.out_filtered_initial_pos_BAM, "w", template=bamfile_2)
    print('FINISHED loading bam files...')
    
    read1_dict = defaultdict(int)
    read2_dict = defaultdict(int)
    bam1_record_dict = defaultdict()
    bam2_record_dict = defaultdict()
    for read1 in bamfile_1:
        read1_record = "{}_{}".format(read1.query_name,read1.query_sequence)
        read1_dict[read1_record] = int(read1.reference_start)
        bam1_record_dict[read1_record] = read1
    
    for read2 in bamfile_2:
        read2_record = "{}_{}".format(read2.query_name,read2.query_sequence)
        read2_dict[read2_record] = int(read2.reference_start)
        bam2_record_dict[read2_record] = read2
        
    bam_intersection_list = list(set(read1_dict.keys()).intersection(set(read2_dict.keys())))
    for read in bam_intersection_list:
        if read1_dict[read] != read2_dict[read]:
            bamfile_diff_filtered.write(bam1_record_dict[read])
            bamfile_inital_pos_filtered.write(bam2_record_dict[read])
            
if __name__ == "__main__":
    sys.exit(main(sys.argv))


