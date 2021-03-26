import matplotlib
matplotlib.use('Agg')
import vcf, argparse, sys
import numpy as np
import pandas as pd
import math
from scipy.stats import chisquare
from collections import defaultdict
import matplotlib.pyplot
import ast

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input 1st resource vcf, 2nd resource vcf, sites vcf containing.')
    parser.add_argument('-a', '--inVCFA', type=str,
        help='Input vcf filepath for sample set A.')
    parser.add_argument('-b', '--inVCFB', type=str,
        help='Input vcf filepath for sample set B.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')

    options = parser.parse_args()
    return options


def main(args):

    options = parse_args()
    
    vcf_reader_A = vcf.Reader(open(options.inVCFA, 'r'))
    vcf_reader_B = vcf.Reader(open(options.inVCFB, 'r'))
    
    
    x_list = list()
    y_list = list()
    for records in vcf.utils.walk_together(vcf_reader_A,vcf_reader_B):
        if None in records: continue
        num_hom_alts = len(record(0).get_hom_alts())
        num_hom_refs = len(record(0).get_hom_refs())
        num_hets = len(record(0).get_hets())
        total_genotypes = float(num_hom_alts + num_hom_refs + num_hets)
        minor_allele_count = (num_hom_alts*2) + num_hets
        minor_af_a = float(minor_allele_count)/float((total_genotypes*2))
        num_hom_alts = len(record(1).get_hom_alts())
        num_hom_refs = len(record(1).get_hom_refs())
        num_hets = len(record(1).get_hets())
        total_genotypes = float(num_hom_alts + num_hom_refs + num_hets)
        minor_allele_count = (num_hom_alts*2) + num_hets
        minor_af_b = float(minor_allele_count)/float((total_genotypes*2))
        x_list.append(minor_af_a)
        y_list.append(minor_af_b)
        
    fig1, ax1 = matplotlib.pyplot.subplots()
    ax1.scatter(x_list, y_list, s=2, color=c_list)
    ax1.set_xlabel('Allele Frequency Sample Set A')
    ax1.set_ylabel('Allele Frequency Sample Set B')
    matplotlib.pyplot.xticks([0.0,0.2,0.4,0.6,0.8,1.0])
    fig1.savefig("AF_comparison.raw_AF.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig1)
   
        
if __name__ == "__main__":
    sys.exit(main(sys.argv))

