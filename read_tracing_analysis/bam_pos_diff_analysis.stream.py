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
    bamfile_1 = pysam.AlignmentFile(options.in_bamFile1, "rb")
    bamfile_2 = pysam.AlignmentFile(options.in_bamFile2, "rb")
    bamfile_diff_filtered = pysam.AlignmentFile(options.out_filtered_moved_pos_BAM, "wb", template=bamfile_1)
    bamfile_inital_pos_filtered = pysam.AlignmentFile(options.out_filtered_initial_pos_BAM, "wb", template=bamfile_2)

    for read1 in bamfile_1:
        read1_name = ""
        if "_1" in read1.query_name.strip() or "/1" in read1.query_name.strip() or read1.is_read1:
            read1_name = "{}/1".format(read1.query_name.strip().strip('_1').strip('/1'))
        elif "_2" in read1.query_name.strip() or "/2" in read1.query_name.strip() or read1.is_read2:
            read1_name = "{}/2".format(read1.query_name.strip().strip('_2').strip('/2'))
        read2_name = "!"
        read1_chr_position = (read1.reference_id, read1.reference_start)
        while read2_name != read1_name:
            if read1_name > read2_name:
                read2 = bamfile_2.next()
            if "_1" in read2.query_name.strip() or "/1" in read2.query_name.strip() or read2.is_read1:
                read2_name = "{}/1".format(read2.query_name.strip().strip('_1').strip('/1'))
            elif "_2" in read2.query_name.strip() or "/2" in read2.query_name.strip() or read2.is_read2:
                read2_name = "{}/2".format(read2.query_name.strip().strip('_2').strip('/2'))
            print('read1_name: {}'.format(read1_name))
            print('read2_name: {}'.format(read2_name))
            
            if read2_name == read1_name:
                read2_chr_position = (read2.reference_id, read2.reference_start)
                if read2_chr_position != read1_chr_position:
                    if (read1_chr_position[0] != read2_chr_position[0]) or abs(int(read1_chr_position[1])-int(read2_chr_position[1])) >= options.pos_diff_thresh:
                        print("{}\t{}:{}\t{}\t{}:{}\t{}\n".format(read1_name,read1_chr_position[0],read1_chr_position[1],read1.mapping_quality,read2_chr_position[0],read2_chr_position[1],read2.mapping_quality), file=options.out_read_pos_diff)
                        if (int(read1.mapping_quality) - int(read2.mapping_quality)) >= options.mapq_diff_thresh:
                            bamfile_diff_filtered.write(read1)
                            bamfile_inital_pos_filtered.write(read2)
                    break
     
if __name__ == "__main__":
    sys.exit(main(sys.argv))


