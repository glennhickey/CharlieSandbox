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
    
    vus_list = set()
    with open(options.inHOMALT, 'r') as homalt_file:
        for line in homalt_file:
            if 'VUS_1|1' in line.split('\t')[1] or 'VUS_1|0' in line.split('\t')[1] or 'VUS_0|1' in line.split('\t')[1]:
                for vus_variant in ast.literal_eval(line.split('\t')[2]):
                    vus_list.update([vus_variant])
    
    vus_variant_genotype_proportions = defaultdict(list)
    non_vus_variant_genotype_proportions = defaultdict(list)
    all_variant_genotype_proportions = defaultdict(list)
    for record in vcf_reader:
        num_hom_alts = len(record.get_hom_alts())
        num_hom_refs = len(record.get_hom_refs())
        num_hets = len(record.get_hets())
        total_genotypes = float(num_hom_alts + num_hom_refs + num_hets)
        hom_ref_frequency = float(num_hom_refs)/float(total_genotypes)
        hom_alt_frequency = float(num_hom_alts)/float(total_genotypes)
        het_frequency = float(num_hets)/float(total_genotypes)
        allele_freq = 1.0 - float(record.INFO['AF'][0])
        variant_record = "{}_{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT[0],allele_freq)
        hom_ref_point = (allele_freq,hom_ref_frequency)
        hom_alt_point = (allele_freq,hom_alt_frequency)
        het_point = (allele_freq,het_frequency)
        all_variant_genotype_proportions['hom_ref'].append(hom_ref_point)
        all_variant_genotype_proportions['hom_alt'].append(hom_alt_point)
        all_variant_genotype_proportions['het'].append(het_point)
        if variant_record in vus_list:
            vus_variant_genotype_proportions['hom_ref'].append(hom_ref_point)
            vus_variant_genotype_proportions['hom_alt'].append(hom_alt_point)
            vus_variant_genotype_proportions['het'].append(het_point)
        else:
            non_vus_variant_genotype_proportions['hom_ref'].append(hom_ref_point)
            non_vus_variant_genotype_proportions['hom_alt'].append(hom_alt_point)
            non_vus_variant_genotype_proportions['het'].append(het_point)
            
     
    x_list = list()
    y_list = list()
    c_list = list()
    for hom_ref_point in vus_variant_genotype_proportions['hom_ref']:
        x_list.append(hom_ref_point(0))
        y_list.append(hom_ref_point(1))
        c_list.append('r')
    for hom_alt_point in vus_variant_genotype_proportions['hom_alt']:
        x_list.append(hom_alt_point(0))
        y_list.append(hom_alt_point(1))
        c_list.append('b')
    for het_point in vus_variant_genotype_proportions['het']:
        x_list.append(het_point(0))
        y_list.append(het_point(1))
        c_list.append('g')
    fig1, ax1 = matplotlib.pyplot.subplots()
    ax1.scatter(x_list, y_list, s=2, color=c_list)
    ax1.set_xlabel('Allele Frequency')
    ax1.set_ylabel('Density')
    fig1.savefig("vus_HWE_AF_distributions.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig1)
    x_list = list()
    y_list = list()
    c_list = list()
    for hom_ref_point in non_vus_variant_genotype_proportions['hom_ref']:
        x_list.append(hom_ref_point(0))
        y_list.append(hom_ref_point(1))
        c_list.append('r')
    for hom_alt_point in non_vus_variant_genotype_proportions['hom_alt']:
        x_list.append(hom_alt_point(0))
        y_list.append(hom_alt_point(1))
        c_list.append('b')
    for het_point in non_vus_variant_genotype_proportions['het']:
        x_list.append(het_point(0))
        y_list.append(het_point(1))
        c_list.append('g')
    fig2, ax2 = matplotlib.pyplot.subplots()
    ax2.scatter(x_list, y_list, s=2, color=c_list)
    ax2.set_xlabel('Allele Frequency')
    ax2.set_ylabel('Density')
    fig2.savefig("non_vus_HWE_AF_distributions.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig2)
    x_list = list()
    y_list = list()
    c_list = list()
    for hom_ref_point in all_variant_genotype_proportions['hom_ref']:
        x_list.append(hom_ref_point(0))
        y_list.append(hom_ref_point(1))
        c_list.append('r')
    for hom_alt_point in all_variant_genotype_proportions['hom_alt']:
        x_list.append(hom_alt_point(0))
        y_list.append(hom_alt_point(1))
        c_list.append('b')
    for het_point in all_variant_genotype_proportions['het']:
        x_list.append(het_point(0))
        y_list.append(het_point(1))
        c_list.append('g')
    fig3, ax3 = matplotlib.pyplot.subplots()
    ax3.scatter(x_list, y_list, s=2, color=c_list)
    ax3.set_xlabel('Allele Frequency')
    ax3.set_ylabel('Density')
    fig3.savefig("all_HWE_AF_distributions.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig3)
     
if __name__ == "__main__":
    sys.exit(main(sys.argv))

