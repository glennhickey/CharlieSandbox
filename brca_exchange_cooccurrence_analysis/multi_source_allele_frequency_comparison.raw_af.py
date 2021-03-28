import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import vcf, argparse, sys
import vcf.utils as vcfutils
import numpy as np
import pandas as pd
import math
from scipy.stats import chisquare
from collections import defaultdict
import matplotlib.pyplot
import ast
import pdb

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input 1st resource vcf, 2nd resource vcf.')
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
    for records in vcfutils.walk_together(vcf_reader_A,vcf_reader_B):
        if None in records:
            print(records)
            continue
        num_hom_alts = len(records[0].get_hom_alts())
        num_hom_refs = len(records[0].get_hom_refs())
        num_hets = len(records[0].get_hets())
        total_genotypes = float(num_hom_alts + num_hom_refs + num_hets)
        minor_allele_count = (num_hom_alts*2) + num_hets
        minor_af_a = float(minor_allele_count)/float((total_genotypes*2))
        num_hom_alts = len(records[1].get_hom_alts())
        num_hom_refs = len(records[1].get_hom_refs())
        num_hets = len(records[1].get_hets())
        total_genotypes = float(num_hom_alts + num_hom_refs + num_hets)
        minor_allele_count = (num_hom_alts*2) + num_hets
        minor_af_b = float(minor_allele_count)/float((total_genotypes*2))
        if minor_af_a != 0.0 or minor_af_b != 0.0:
            if minor_af_a == 0.0: minor_af_a = 1.0
            if minor_af_b == 0.0: minor_af_b = 1.0
            x_list.append(minor_af_a)
            y_list.append(minor_af_b)
    
    adjust_figure = False
    if min(x_list) < 0.001:
        print("WARNING an x value is less than 0.001")
        adjust_figure = True
    
    if min(y_list) < 0.001:
        print("WARNING an y value is less than 0.001")
        adjust_figure = True
    
    if adjust_figure: pdb.set_trace()
    fig = plt.figure()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.scatter(x_list, y_list, color='black', s=8)
    ax.set_xlabel('Allele Frequency Sample Set Broad')
    ax.set_ylabel('Allele Frequency Sample Set Not Broad')
    ax.grid()
    ax.plot([0.0, 1.0], [0.0, 1.0], color="red")
    plt.xticks([0.001,0.01,0.1,1])
    plt.yticks([0.001,0.01,0.1,1])
    plt.title("Allele frequency for {}".format(options.outReport))
    fig.savefig("AF_comparison.raw_AF.{}.png".format(options.outReport))
    plt.close(fig)
        
if __name__ == "__main__":
    sys.exit(main(sys.argv))

