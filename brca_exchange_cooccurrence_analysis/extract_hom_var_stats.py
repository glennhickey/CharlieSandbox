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
    parser = argparse.ArgumentParser('Input VCF and output homozygous alt report and histogram.')
    parser.add_argument('-i', '--inVCF', type=str,
        help='Input vcf filepath.')
    parser.add_argument('-v', '--inHOMALT', type=str,
        help='Input list of samples that contain at least 1 homozygous alt VUS genotype.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')

    options = parser.parse_args()
    return options


def main(args):

    options = parse_args()
    
    vcf_reader = vcf.Reader(open(options.inVCF, 'r'))
    hom_var_dict = defaultdict(int)
    het_dict = defaultdict(int)
    hom_ref_dict = defaultdict(int)
    hom_var_af_dict = defaultdict(list)
    het_af_dict = defaultdict(list)
    hom_ref_af_dict = defaultdict(list)
    
    for record in vcf_reader:
        for call in record.get_hom_alts():
            if not hom_var_dict[call.sample]: hom_var_dict[call.sample] = 0
            hom_var_dict[call.sample] += 1
            hom_var_af_dict[call.sample].append(record.INFO['AF'][0])
        for call in record.get_hets():
            if not het_dict[call.sample]: het_dict[call.sample] = 0
            het_dict[call.sample] += 1
            het_af_dict[call.sample].append(record.INFO['AF'][0])
        for call in record.get_hom_refs():
            if not hom_ref_dict[call.sample]: hom_ref_dict[call.sample] = 0
            hom_ref_dict[call.sample] += 1
            hom_ref_af_dict[call.sample].append(record.INFO['AF'][0])
    
    hom_vus_list = defaultdict(int)
    hom_vus_af_list = defaultdict(list)
    with open(options.inHOMALT, 'r') as homalt_file:
        for line in homalt_file:
            if 'VUS_1|1' in line.split('\t')[1]:
                sample_id = line.split('\t')[0]
                variant_list = ast.literal_eval(line.split('\t')[2])
                hom_vus_list[sample_id] = len(variant_list)
                for variant in variant_list:
                    variant_af = variant.split('_')[4]
                    hom_vus_af_list[sample_id].append(variant_af)
    
    import pdb; pdb.set_trace()
    
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
    fig1, ax1 = matplotlib.pyplot.subplots()
    logbins = np.geomspace(min_genotypes, max_genotypes, 100)
    ax1.hist(hom_var_dict.values(), bins=logbins)
    ax1.set_title("Hom ALT Genotype Sample Distribution LOG")
    matplotlib.pyplot.xscale('log')
    fig1.savefig("hom_alt_dist_log.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig1)
    
    fig2, ax2 = matplotlib.pyplot.subplots()
    ax2.hist(hom_var_dict.values(), bins=100)
    ax2.set_title("Hom ALT Genotype Sample Distribution LINEAR")
    matplotlib.pyplot.xscale('linear')
    fig2.savefig("hom_alt_dist_linear.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig2)
    
    # hom vs het proportion analysis
    hom_prop_dict = defaultdict(float)
    het_prop_dict = defaultdict(float)
    for sample in hom_var_dict.keys():
        hom_prop_dict[sample] = float(hom_var_dict[sample]) / (float(hom_var_dict[sample]) + float(het_dict[sample]) + float(hom_ref_dict[sample]))
        het_prop_dict[sample] = float(het_dict[sample]) / (float(hom_var_dict[sample]) + float(het_dict[sample]) + float(hom_ref_dict[sample]))
    
    x_list = list()
    y_list = list()
    c_list = list()
    for sample in hom_var_dict.keys():
        if sample in hom_vus_list:
            x_list.append(hom_vus_list[sample])
            y_list.append(het_dict[sample])
            c_list.append('r')
    fig3, ax3 = matplotlib.pyplot.subplots()
    ax3.scatter(x_list, y_list, s=2, color=c_list)
    ax3.set_title("Homozygous VUS Counts to Heterozygous Counts")
    ax3.set_xlabel('homozygous VUS counts')
    ax3.set_ylabel('heterozygous counts')
    fig3.savefig("hom_VUS_het_counts_scatter.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig3)
    x_list = list()
    y_list = list()
    c_list = list()
    for sample in hom_var_dict.keys():
        if sample in hom_vus_list:
            x_list.append(hom_var_dict[sample] - hom_vus_list[sample])
            y_list.append(het_dict[sample])
            c_list.append('b')
    fig4, ax4 = matplotlib.pyplot.subplots()
    ax4.scatter(x_list, y_list, s=2, color=c_list)
    ax4.set_title("Homozygous non-VUS Counts to Heterozygous Counts")
    ax4.set_xlabel('homozygous non-VUS counts')
    ax4.set_ylabel('heterozygous counts')
    fig4.savefig("hom_non_VUS_het_counts_scatter.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig4)
     
if __name__ == "__main__":
    sys.exit(main(sys.argv))

