import vcf, argparse, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from scipy.stats import chisquare
from collections import defaultdict

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input VCF and output homozygous alt report and histogram.')
    parser.add_argument('-i', '--inVCF', type=str,
        help='Input vcf filepath.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')

    options = parser.parse_args()
    return options


def main(args):

    options = parse_args()
    
    vcf_reader = vcf.Reader(open(options.inVCF, 'r'))
    hom_var_dict = defaultdict(int)

    for record in vcf_reader:
        for call in record.get_hom_alts():
            hom_var_dict[call.sample] += 1
    
    # Compile and output report
    num_samples = len(hom_var_dict)
    min_genotypes = min(hom_var_dict.values())
    Q1_genotypes = list(sorted(hom_var_dict.values()))[-int(len(hom_var_dict.values())*0.75)]
    median_genotypes = list(sorted(hom_var_dict.values()))[-int(len(hom_var_dict.values())*0.5)]
    Q3_genotypes = list(sorted(hom_var_dict.values()))[-int(len(hom_var_dict.values())*0.25)]
    max_genotypes = max(hom_var_dict.values())
    with open(options.outReport, 'w') as report_file:
        report_file.write("Number of samples: {}\n".format(num_samples))
        report_file.write("Number of homozygous alt genotypes for a sample:\n")
        report_file.write("Min: {}\n".format(min_genotypes))
        report_file.write("1st Quartile: {}\n".format(Q1_genotypes))
        report_file.write("Median: {}\n".format(median_genotypes))
        report_file.write("3rd Quartile: {}\n".format(Q3_genotypes))
        report_file.write("Max: {}\n".format(max_genotypes))
        
    
    # Build histogram plots and save to .png files
    fig1, ax1 = plt.subplots()
    logbins = np.geomspace(min_genotypes, max_genotypes, 100)
    ax1.hist(hom_var_dict.values(), bins=logbins)
    ax1.set_title("Hom ALT Genotype Sample Distribution LOG")
    plt.xscale('log')
    fig1.savefig("hom_alt_dist_log.{}.png".format(options.outReport))
    plt.close(fig1)
    
    fig2, ax2 = plt.subplots()
    ax2.hist(hom_var_dict.values(), bins=100)
    ax2.set_title("Hom ALT Genotype Sample Distribution LINEAR")
    plt.xscale('linear')
    fig2.savefig("hom_alt_dist_linear.{}.png".format(options.outReport))
    plt.close(fig2)
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))

