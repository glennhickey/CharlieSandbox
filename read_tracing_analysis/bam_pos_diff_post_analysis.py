#!/usr/bin/env python2.7

"""
Charles Markello
2/12/2019

Description:
    

"""

from __future__ import print_function
import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_args():
    """
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """

    parser = argparse.ArgumentParser(prog='read-tracer-mapq-plots', description= __doc__)
    parser.add_argument('in_sorted_bam_diff_file', nargs='?', type=str)
    
    options = parser.parse_args()

    return options

def main(args):
    
    options = parse_args()

    input_file = open(options.in_sorted_bam_diff_file, 'r')
    count_total_records = 0
    count_nomap_nomap = 0
    count_nomap_map = 0
    nomap_map_mapq_list = list()
    count_map_nomap = 0
    map_nomap_mapq_list = list()
    count_map_map = 0
    map_map_pos_diff_list = list()
    map_map_mapq_1st_greater_list = list()
    map_map_mapq_2nd_greater_list = list()
    map_map_mapq_equal_list = list()
    map_map_mapq_combined_list_total = list()
    map_map_mapq_combined_list_no_mapq0 = list() 

    for record_diff in input_file:
        count_total_records = count_total_records + 1
        record_diff_list = record_diff.strip().split('\t')
        
        if record_diff_list[1] == '-1' and record_diff_list[4] == '-1':
            count_nomap_nomap = count_nomap_nomap + 1
        elif record_diff_list[1] == '-1' and record_diff_list[4] != '-1':
            count_nomap_map = count_nomap_map + 1
            nomap_map_mapq_list.append(int(record_diff_list[6]))
        elif record_diff_list[1] != '-1' and record_diff_list[4] == '-1':
            count_map_nomap = count_map_nomap + 1
            map_nomap_mapq_list.append(int(record_diff_list[3]))
        elif record_diff_list[1] != '-1' and record_diff_list[4] != '-1' and record_diff_list[1] == record_diff_list[4]:
            count_map_map = count_map_map + 1
            bam1_position = int(record_diff_list[2])
            bam1_mapq = int(record_diff_list[3])
            bam2_position = int(record_diff_list[5])
            bam2_mapq = int(record_diff_list[6])
            map_position_difference = abs(bam1_position - bam2_position)
            map_map_pos_diff_list.append(map_position_difference)
            
            map_map_mapq_combined_list_total.append(bam1_mapq - bam2_mapq)
            
            if bam1_mapq > bam2_mapq:
                map_map_mapq_1st_greater_list.append(bam1_mapq - bam2_mapq)
            elif bam1_mapq < bam2_mapq:
                map_map_mapq_2nd_greater_list.append(bam2_mapq - bam1_mapq)
            elif bam1_mapq == bam2_mapq:
                map_map_mapq_equal_list.append(bam1_mapq)
            
            if bam1_mapq != bam2_mapq:
                map_map_mapq_combined_list_no_mapq0.append(bam1_mapq - bam2_mapq)    
    
    hist_map_map_mapq_1st_greater = np.histogram(map_map_mapq_1st_greater_list, range=(1,60), bins=59)
    hist_map_map_mapq_2nd_greater = np.histogram(map_map_mapq_2nd_greater_list, range=(1,60), bins=59)
    fig = plt.figure()
    plot_figure = plt.hist(map_map_mapq_1st_greater_list, range=(1,60), bins=59)
    fig.savefig('hist_map_map_mapq_1st_greater.png')
    fig = plt.figure()
    plot_figure = plt.hist(map_map_mapq_2nd_greater_list, range=(1,60), bins=59)
    fig.savefig('hist_map_map_mapq_2nd_greater.png')
    fig = plt.figure()
    #plot_figure = plt.hist(map_map_mapq_1st_greater_list, map_map_mapq_2nd_greater_list, range=(1,60), bins=[59, 59])
    #fig.savefig('hist_map_map_mapq_both.png')
    plot_figure = plt.hist(map_map_mapq_combined_list_total, range=(-60,60), bins=120)
    plt.xlabel('mapq difference')
    fig.savefig('hist_map_map_mapq_combined_total.png')
    plot_figure = plt.hist(map_map_mapq_combined_list_no_mapq0, range=(-60,60), bins=119)
    plt.xlabel('mapq difference')
    fig.savefig('hist_map_map_mapq_combined_no_mapq0.png') 
 
if __name__ == "__main__":
    sys.exit(main(sys.argv))


