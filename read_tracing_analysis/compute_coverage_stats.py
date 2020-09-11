#!/usr/bin/env python2.7

"""
Charles Markello

Description:
    Extracts coverage information from input bam files, compares them and generates plots displaying relative coverage statistics

Assumptions:
    All exon positions for a given gene that's listed in the input bed file lies on the same contig.    

"""

from __future__ import print_function
import argparse
import sys
from collections import defaultdict
import pysam

def parse_args():
    """
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """

    parser = argparse.ArgumentParser(prog='read-tracer', description= __doc__)
    parser.add_argument('-i', '--in_bam_file_1', type=str,
        help="Input 1st bam file")
    parser.add_argument('-j', '--in_bam_file_2', type=str,
        help="Input 2nd bam file")
    parser.add_argument('-b', '--in_gene_bed_file', type=str,
        help="Input gene bed file. 4th column must have gene name.")
    parser.add_argument('-o','--out_file_basename', type=str, 
        help="Basename for output files.")
    
    options = parser.parse_args()

    return options

def main(args):
    
    options = parse_args()
    
    # Extract gene boundaries from exon-level gene bed files
    gene_exon_boundary_dict = defaultdict(list)
    with open(options.in_gene_bed_file, 'r') as gene_bed_file:
        for line in gene_bed_file:
            split_line = line.strip().split()
            chromosome = split_line[0]
            start_pos = split_line[1]
            end_pos = split_line[2]
            gene_name = split_line[3]
            gene_exon_boundary_dict[gene_name].append('{}:{}-{}'.format(chromosome,start_pos,end_pos))
    
    gene_boundary_dict = defaultdict(str)
    for gene_name in gene_exon_boundary_dict.keys():
        start_pos_list = list()
        end_pos_list = list()
        chromosome_pos = ""
        for exon_boundary in gene_exon_boundary_dict[gene_name]:
            split_exon_boundary = exon_boundary.split(':')
            chromosome_pos = split_exon_boundary[0]
            start_pos_list.append(int(split_exon_boundary[1].split('-')[0]))
            end_pos_list.append(int(split_exon_boundary[1].split('-')[1]))
        min_start_pos = min(start_pos_list)
        max_end_pos = max(end_pos_list)
        gene_boundary_dict[gene_name] = '{}:{}-{}'.format(chromosome_pos,min_start_pos,max_end_pos)
    
    # Extract exon-level coverage stats in bam files
    bam_exon_coverage_dict = defaultdict(tuple)
    for gene_name in gene_exon_boundary_dict.keys():
        exon_boundary_list = gene_exon_boundary_dict[gene_name]
        exon_id = 0
        for exon_boundary in exon_boundary_list:
            bam_1_coverage_stats = pysam.coverage("-r", exon_boundary, options.in_bam_file_1).strip().split('\n')[1].split('\t')
            bam_2_coverage_stats = pysam.coverage("-r", exon_boundary, options.in_bam_file_2).strip().split('\n')[1].split('\t')
            print(bam_1_coverage_stats)
            print(bam_2_coverage_stats)
            print()
            total_region_length = float(int(bam_1_coverage_stats[2]) - int(bam_1_coverage_stats[1])+1)
            bam_1_covered_bases = float(bam_1_coverage_stats[4])
            bam_1_percent_covered = float(bam_1_covered_bases/total_region_length)*100.00
            bam_2_covered_bases = float(bam_2_coverage_stats[4])
            bam_2_percent_covered = float(bam_2_covered_bases/total_region_length)*100.00
            bam_exon_coverage_dict['{}_{}'.format(gene_name,exon_id)] = (bam_1_percent_covered,bam_2_percent_covered)
            exon_id += 1
    
    for exon_name in bam_exon_coverage_dict.keys():
        print('exon_name: {}. Relative coverage percentage {}'.format(exon_name,bam_exon_coverage_dict[exon_name]))
        
    # Extract gene-level coverage stats in bam files
    bam_gene_coverage_dict = defaultdict(tuple)
    for gene_name in gene_boundary_dict.keys():
        gene_boundary = gene_boundary_dict[gene_name]
        bam_1_coverage_stats = pysam.coverage("-r", gene_boundary, options.in_bam_file_1).strip().split('\n')[1].split('\t')
        bam_2_coverage_stats = pysam.coverage("-r", gene_boundary, options.in_bam_file_2).strip().split('\n')[1].split('\t')
        print(bam_1_coverage_stats)
        print(bam_2_coverage_stats)
        total_region_length = float(int(bam_1_coverage_stats[2]) - int(bam_1_coverage_stats[1])+1)
        bam_1_covered_bases = float(bam_1_coverage_stats[4])
        bam_1_percent_covered = float(bam_1_covered_bases/total_region_length)*100.00
        bam_2_covered_bases = float(bam_2_coverage_stats[4])
        bam_2_percent_covered = float(bam_2_covered_bases/total_region_length)*100.00
        bam_gene_coverage_dict[gene_name] = (bam_1_percent_covered,bam_2_percent_covered)
    
    for gene_name in bam_gene_coverage_dict.keys():
        print('gene_name: {}. Relative coverage percentage {}'.format(gene_name,bam_gene_coverage_dict[gene_name]))
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))


