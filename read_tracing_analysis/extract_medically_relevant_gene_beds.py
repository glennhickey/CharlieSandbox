#!/usr/bin/env python2.7

"""
Charles Markello

Description:
    

"""

from __future__ import print_function
import argparse
import sys
from collections import defaultdict
import csv

def parse_args():
    """
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """

    parser = argparse.ArgumentParser(prog='read-tracer', description= __doc__)
    parser.add_argument('-n', '--in_gene_csv', type=str,
        help="Input NGS-level gene list CSV file")
    parser.add_argument('-e', '--in_exon_csv', type=str,
        help="Input NGS-level exon list CSV file")
    parser.add_argument('-a', '--in_association_level', type=int,
        help="Input minimum gene-disease association level (3,2,1,0). 3=definite, 2=likely, 1=weakly, 0=undetermined.")
    parser.add_argument('-o','--out_bed_file', type=str, 
        help="Output bed file of gene regions.")
    
    options = parser.parse_args()

    return options

def main(args):
    
    options = parse_args()
    
    gene_list_association_level_list = list()
    with open(options.in_gene_csv, 'r') as gene_list_csv:
        gene_list_csv_reader = csv.reader(gene_list_csv, delimiter=',', quotechar='"')
        for line in gene_list_csv_reader:
            gene_name = line[0]
            if gene_name == 'gene': continue
            print(line)
            if line[13] == '':
                association_level = 0
            else:
                association_level = int(line[13])
            if association_level is None: association_level = 0
            if association_level >= options.in_association_level:
                gene_list_association_level_list.append(gene_name)
    
    print(len(gene_list_association_level_list))
    with open(options.in_exon_csv, 'r') as exon_list_csv, open(options.out_bed_file, 'w') as output_bed_file:
        exon_list_csv_reader = csv.reader(exon_list_csv, delimiter=',', quotechar='"')
        for line in exon_list_csv_reader:
            chromosome = line[0]
            if chromosome == 'chr': continue
            start_minus_65bp = int(line[1]) - 1
            end_plit_65bp = int(line[2])
            gene_name = line[5]
            if gene_name in gene_list_association_level_list:
                output_bed_file.write('{}\t{}\t{}\t{}\n'.format(chromosome, start_minus_65bp, end_plit_65bp, gene_name))
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))


