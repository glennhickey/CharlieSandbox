import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import vcf, argparse, sys
import vcf.utils as vcfutils
import numpy as np
import math
from scipy.stats import ks_2samp
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
    rare_x_list = list()
    rare_y_list = list()
    common_x_list = list()
    common_y_list = list()
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
        if minor_af_a == 0.0 or minor_af_b == 0.0: continue
        x_list.append(minor_af_a)
        y_list.append(minor_af_b)
        if minor_af_a == 0.01 or minor_af_b <= 0.01:
            rare_x_list.append(minor_af_a)
            rare_y_list.append(minor_af_b)
        else:
            common_x_list.append(minor_af_a)
            common_y_list.append(minor_af_b)
        #if minor_af_a != 0.0 or minor_af_b != 0.0:
        #    if minor_af_a == 0.0: minor_af_a = 1.0
        #    if minor_af_b == 0.0: minor_af_b = 1.0
        #    x_list.append(minor_af_a)
        #    y_list.append(minor_af_b)
    
    adjust_figure = False
    if min(x_list) < 0.001:
        print("WARNING an x value is less than 0.001")
        adjust_figure = True
    
    if min(y_list) < 0.001:
        print("WARNING an y value is less than 0.001")
        adjust_figure = True
    
    # Calculate ks stats
    total_ks_stat,total_ks_pvalue = ks_2samp(x_list, y_list)
    rare_ks_stat,rare_ks_pvalue = ks_2samp(rare_x_list, rare_y_list)
    common_ks_stat,common_ks_pvalue = ks_2samp(common_x_list, common_y_list)
    with open('AF_comparison.raw_AF.{}.ks_stats.tsv'.format(options.outReport), 'w') as ks_output:
        ks_output.write('variant_type\tsample_size\tks_stat\tks_pvalue\n')
        ks_output.write('total\t{}\t{}\t{}\n'.format(len(x_list),total_ks_stat,total_ks_pvalue))
        ks_output.write('rare\t{}\t{}\t{}\n'.format(len(rare_x_list),rare_ks_stat,rare_ks_pvalue))
        ks_output.write('common\t{}\t{}\t{}\n'.format(len(common_x_list),common_ks_stat,common_ks_pvalue))
    
    if adjust_figure: pdb.set_trace()
    # Plot total variant figure
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
    
    # Plot rare variant figure
    fig = plt.figure()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.scatter(rare_x_list, rare_y_list, color='black', s=8)
    ax.set_xlabel('Allele Frequency Sample Set Broad')
    ax.set_ylabel('Allele Frequency Sample Set Not Broad')
    ax.grid()
    ax.plot([0.0, 1.0], [0.0, 1.0], color="red")
    plt.xticks([0.001,0.01,0.1,1])
    plt.yticks([0.001,0.01,0.1,1])
    plt.title("Allele frequency for {}".format(options.outReport))
    fig.savefig("AF_comparison.raw_AF.either_less_than_or_equal_0.01.{}.png".format(options.outReport))
    plt.close(fig)
    
    # Plot common variant figure
    fig = plt.figure()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.scatter(common_x_list, common_y_list, color='black', s=8)
    ax.set_xlabel('Allele Frequency Sample Set Broad')
    ax.set_ylabel('Allele Frequency Sample Set Not Broad')
    ax.grid()
    ax.plot([0.0, 1.0], [0.0, 1.0], color="red")
    plt.xticks([0.001,0.01,0.1,1])
    plt.yticks([0.001,0.01,0.1,1])
    plt.title("Allele frequency for {}".format(options.outReport))
    fig.savefig("AF_comparison.raw_AF.both_more_than_0.01.{}.png".format(options.outReport))
    plt.close(fig)
        
if __name__ == "__main__":
    sys.exit(main(sys.argv))

