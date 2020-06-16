import matplotlib
matplotlib.use('Agg')
import vcf, argparse, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from scipy.stats import chisquare
from collections import defaultdict
import hail as hl

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input both a VUS and pathogenic concordant bgzip compressed and tabix indexed vcf.')
    parser.add_argument('-i', '--inVUSvcf', type=str,
        help='Input brcaexchange-VUS/sample-genotype concordant vcf filepath.')
    parser.add_argument('-j', '--inPATHvcf', type=str,
        help='Input brcaexchange-pathogenic/sample-genotype concordant vcf filepath.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')
    parser.add_argument('-v', '--outVariants', type=str,
        help='Output apparent-benign VUS variants list filename.')

    options = parser.parse_args()
    return options

# Shamlessly borrowed from http://specminor.org/2017/01/08/performing-chi-squared-gof-python.html
# Becuase I need to calculate chi-square p-values at greater resolution than what's provided
# by scipy.stats.chisquare
def gf(x):
    #Play with these values to adjust the error of the approximation.
    upper_bound=100.0
    resolution=1000000.0

    step=upper_bound/resolution

    val=0
    rolling_sum=0

    while val<=upper_bound:
        val+=step
        rolling_sum+=step*(val**(x-1)*2.7182818284590452353602874713526624977**(-val))

    return rolling_sum

def ilgf(s,z):
    val=0

    for k in range(0,100):
        val+=(((-1)**k)*z**(s+k))/(math.factorial(k)*(s+k))
    return val

def chisquarecdf(x,k):
    return 1-ilgf(k/2,x/2)/gf(k/2)

def chisquare_custom(observed_values,expected_values):
    test_statistic=0

    for observed, expected in zip(observed_values, expected_values):
        test_statistic+=(abs(float(observed)-float(expected))-0.5)**2/float(expected)

    df=1
    print("test_statistic: {}".format(test_statistic))
    return test_statistic, chisquarecdf(test_statistic,df)

def main(args):

    options = parse_args()

    ## Isolate samples by the 2 different phased PATHOGENIC het calls
    with open(options.inPATHvcf, 'rb') as inPATHvcf_file:
        vcf_reader_pathogenic = vcf.Reader(inPATHvcf_file)
        brca_pathogenic_left_het_sample_list = defaultdict(list)
        brca_pathogenic_right_het_sample_list = defaultdict(list)
        brca_pathogenic_hom_sample_list = defaultdict(list)
        for record in vcf_reader_pathogenic:
            variant_record = "{}_{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT[0],record.INFO['AF'][0])
            for sample in record.samples:
                if sample['GT'] == '1|0':
                    brca_pathogenic_left_het_sample_list[sample.sample].append(variant_record)
                elif sample['GT'] == '0|1':
                    brca_pathogenic_right_het_sample_list[sample.sample].append(variant_record)
                elif sample['GT'] == '1|1':
                    brca_pathogenic_hom_sample_list[sample.sample].append(variant_record)
    
    ## Isolate samples by the 2 different phased VUS het calls
    VUS_hom_HWE_stats = defaultdict(list)
    with open(options.inVUSvcf, 'rb') as inVUSvcf_file:
        vcf_reader_vus = vcf.Reader(inVUSvcf_file)
        brca_vus_left_het_sample_list = defaultdict(list)
        brca_vus_right_het_sample_list = defaultdict(list)
        brca_vus_hom_sample_list = defaultdict(list)
        for record in vcf_reader_vus:
            variant_record = "{}_{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT[0],record.INFO['AF'][0])
            HWE_obs_genotype_freq = list()
            HWE_exp_genotype_freq = list()
            HWE_obs_genotype_freq.append(0)
            HWE_obs_genotype_freq.append(0)
            HWE_obs_genotype_freq.append(0)
            print(variant_record)
            for sample in record.samples:
                if sample['GT'] == '0|0':
                    HWE_obs_genotype_freq[0] += 1
                elif sample['GT'] == '1|0':
                    HWE_obs_genotype_freq[1] += 1
                    brca_vus_left_het_sample_list[sample.sample].append(variant_record)
                elif sample['GT'] == '0|1':
                    HWE_obs_genotype_freq[1] += 1
                    brca_vus_right_het_sample_list[sample.sample].append(variant_record)
                elif sample['GT'] == '1|1':
                    HWE_obs_genotype_freq[2] += 1
                    brca_vus_hom_sample_list[sample.sample].append(variant_record)
            # Hardy Weinberg Calculation
            if HWE_obs_genotype_freq[2] > 0:
                q_allele_freq = float(record.INFO['AF'][0])
                p_allele_freq = 1.0 - float(record.INFO['AF'][0])
                HWE_exp_genotype_freq.append((p_allele_freq * p_allele_freq)*len(record.samples))
                HWE_exp_genotype_freq.append((2.0 * p_allele_freq * q_allele_freq)*len(record.samples))
                HWE_exp_genotype_freq.append((q_allele_freq * q_allele_freq)*len(record.samples))
                chisquare_value = chisquare(HWE_obs_genotype_freq,HWE_exp_genotype_freq,ddof=1)
                levene_haldane_hwe_calc = hl.eval(hl.hardy_weinberg_test(HWE_obs_genotype_freq[0],HWE_obs_genotype_freq[1],HWE_obs_genotype_freq[2]))
                VUS_hom_HWE_stats[variant_record] = [HWE_obs_genotype_freq, HWE_exp_genotype_freq, p_allele_freq, q_allele_freq, chisquare_value[0], chisquare_value[1], levene_haldane_hwe_calc.het_freq_hwe, levene_haldane_hwe_calc.p_value]
    
    ## Output to hom var VUS Hardy-Weinberg Equilibrium report
    hwe_report_filename = "hom_vus_hwe_{}".format(options.outReport)
    hwe_chi_square_stats = list()
    hwe_chi_square_pvalues = list()
    minor_allele_freqs = list()
    hwe_hl_dist_pvalues = list()
    with open(hwe_report_filename, 'w') as hwe_report_file:
        hwe_report_file.write("variant_record\thwe_obs_(0/0,0/1,1/1)\thwe_exp_(0/0,0/1,1/1)\tp_freq\tq_freq\tscipy_chi_square_stat\tscipy_chi_square_p_value\n")   
        for variant_record in VUS_hom_HWE_stats.keys():
            hwe_stats = VUS_hom_HWE_stats[variant_record]
            hwe_chi_square_stats.append(hwe_stats[4])
            hwe_chi_square_pvalues.append(hwe_stats[5])
            minor_allele_freqs.append(hwe_stats[3])
            hwe_hl_dist_pvalues.append(hwe_stats[7])
            hwe_report_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(variant_record, hwe_stats[0], hwe_stats[1], hwe_stats[2], hwe_stats[3], hwe_stats[4], hwe_stats[5], hwe_stats[6], hwe_stats[7]))
    # Build histogram plots and save to .png files
    print("hwe_chi_square_stats: {}".format(hwe_chi_square_stats))
    fig1, ax1 = plt.subplots()
    hwe_chi_square_stats = pd.Series(hwe_chi_square_stats)
    logbins = np.geomspace(hwe_chi_square_stats.min(), hwe_chi_square_stats.max(), 100)
    ax1.hist(hwe_chi_square_stats, bins=logbins)
    ax1.set_title("Hom VUS chi square stat distribution")
    plt.xscale('log')
    fig1.savefig("hom_vus_chi_square_stat.{}.png".format(options.outReport))
    plt.close(fig1)
    print("hwe_chi_square_pvalues: {}".format(hwe_chi_square_pvalues))
    fig2, ax2 = plt.subplots()
    ax2.hist(hwe_chi_square_pvalues, bins=50)
    ax2.set_title("Hom VUS corrected chi square pvalue distribution")
    plt.xscale('linear')
    fig2.savefig("hom_vus_chi_square_pvalue.{}.png".format(options.outReport))
    plt.close(fig2)
    print("minor_allele_freqs: {}".format(minor_allele_freqs))
    fig3, ax3 = plt.subplots()
    minor_allele_freqs = pd.Series(minor_allele_freqs)
    logbins = np.geomspace(minor_allele_freqs.min(), minor_allele_freqs.max(), 100)
    ax3.hist(minor_allele_freqs, bins=logbins)
    ax3.set_title("Hom VUS allele frequency distribution")
    plt.xscale('log')
    fig3.savefig("hom_vus_allele_frequencies.{}.png".format(options.outReport))
    plt.close(fig3)
    print("hwe_hl_dist_pvalues: {}".format(hwe_hl_dist_pvalues))
    fig4, ax4 = plt.subplots()
    ax4.hist(hwe_hl_dist_pvalues, bins=50)
    ax4.set_title("Hom VUS Levene-Haldane exact test pvalue distribution")
    plt.xscale('linear')
    fig4.savefig("hom_vus_hl_exact_pvalue.{}.png".format(options.outReport))
    plt.close(fig4)
    import pdb; pdb.set_trace()

    ## Single genotype samples sets ##
    # samples containing PATH 1|0 genotype
    brca_pathogenic_left_het_sample_set = set(brca_pathogenic_left_het_sample_list.keys())
    
    # samples containing PATH 0|1 genotype
    brca_pathogenic_right_het_sample_set = set(brca_pathogenic_right_het_sample_list.keys())
    
    # samples containing PATH 1|1 genotype
    brca_pathogenic_hom_sample_set = set(brca_pathogenic_hom_sample_list.keys())
    
    # samples containing VUS 1|0 genotype
    brca_vus_left_het_sample_set = set(brca_vus_left_het_sample_list.keys())
    
    # samples containing VUS 0|1 genotype
    brca_vus_right_het_sample_set = set(brca_vus_right_het_sample_list.keys())
    
    # samples containing VUS 1|1 genotype
    brca_vus_hom_sample_set = set(brca_vus_hom_sample_list.keys())

    ## Compound genotype sample sets ##
    # samples containing both PATH 1|0 and PATH 0|1 genotypes
    brca_left_het_path_right_het_path_set = brca_pathogenic_left_het_sample_set.intersection(brca_pathogenic_right_het_sample_set)
    
    # samples containing both PATH 1|0 and VUS 1|0 genotypes
    brca_left_het_coocourance_sample_set = brca_pathogenic_left_het_sample_set.intersection(brca_vus_left_het_sample_set)
    
    # samples containing both PATH 1|0 and VUS 0|1 genotypes
    brca_left_het_path_right_het_vus_coocourance_sample_set = brca_pathogenic_left_het_sample_set.intersection(brca_vus_right_het_sample_set)
    
    # samples containing both PATH 1|0 and VUS 1|1 genotypes
    brca_left_het_path_hom_vus_coocourance_sample_set = brca_pathogenic_left_het_sample_set.intersection(brca_vus_hom_sample_set)
    
    # samples containing both PATH 0|1 and VUS 1|0 genotypes
    brca_right_het_path_left_het_vus_coocourance_sample_set = brca_pathogenic_right_het_sample_set.intersection(brca_vus_left_het_sample_set)
    
    # samples containing both PATH 0|1 and VUS 0|1 genotypes
    brca_right_het_path_right_het_vus_coocourance_sample_set = brca_pathogenic_right_het_sample_set.intersection(brca_vus_right_het_sample_set)
    
    # samples containing both PATH 0|1 and VUS 1|1 genotypes
    brca_right_het_path_hom_vus_coocourance_sample_set = brca_pathogenic_right_het_sample_set.intersection(brca_vus_hom_sample_set)
    
    # samples containing both PATH 1|1 and VUS 1|0 genotypes
    brca_hom_path_left_het_vus_coocourance_sample_set = brca_pathogenic_hom_sample_set.intersection(brca_vus_left_het_sample_set)
    
    # samples containing both PATH 1|1 and VUS 0|1 genotypes
    brca_hom_path_right_het_vus_coocourance_sample_set = brca_pathogenic_hom_sample_set.intersection(brca_vus_right_het_sample_set)
    
    # samples containing both PATH 1|1 and VUS 1|1 genotypes
    brca_hom_path_hom_vus_coocourance_sample_set = brca_pathogenic_hom_sample_set.intersection(brca_vus_hom_sample_set)

    # samples containing either (PATH 1|0 and VUS 0|1), (PATH 0|1 and VUS 1|0) or (VUS 1|1)
    brca_trans_union_sample_set = set.union(brca_left_het_path_right_het_vus_coocourance_sample_set, brca_right_het_path_left_het_vus_coocourance_sample_set, brca_vus_hom_sample_set)

    ## Look at unique vus variants that coincide with concurrance sample set
    # VUS variants in samples with (PATH 1|0 and VUS 0|1)
    brca_concurrent_vus_coordinates_category_8 = list()
    for sample in brca_left_het_path_right_het_vus_coocourance_sample_set:
        brca_concurrent_vus_coordinates_category_8.append(brca_vus_right_het_sample_list[sample])

    # VUS variants in samples with (PATH 1|0 and VUS 1|1)
    brca_concurrent_vus_coordinates_category_9 = list()
    for sample in brca_left_het_path_hom_vus_coocourance_sample_set:
        brca_concurrent_vus_coordinates_category_9.append(brca_vus_hom_sample_list[sample])

    # VUS variants in samples with (PATH 0|1 and VUS 1|0)
    brca_concurrent_vus_coordinates_category_10 = list()
    for sample in brca_right_het_path_left_het_vus_coocourance_sample_set:
        brca_concurrent_vus_coordinates_category_10.append(brca_vus_left_het_sample_list[sample])

    # VUS variants in samples with (PATH 0|1 and VUS 1|1)
    brca_concurrent_vus_coordinates_category_12 = list()
    for sample in brca_right_het_path_hom_vus_coocourance_sample_set:
        brca_concurrent_vus_coordinates_category_12.append(brca_vus_hom_sample_list[sample])
    
    # VUS variants in samples with (VUS 1|1)
    brca_concurrent_vus_coordinates_category_13 = list()
    for sample in brca_vus_hom_sample_set:
        brca_concurrent_vus_coordinates_category_13.append(brca_vus_hom_sample_list[sample])
    
    
    # Flatten out the variant lists into sets
    brca_concurrent_vus_coordinates_category_8_set = set([val for sublist in brca_concurrent_vus_coordinates_category_8 for val in sublist])
    brca_concurrent_vus_coordinates_category_9_set = set([val for sublist in brca_concurrent_vus_coordinates_category_9 for val in sublist])
    brca_concurrent_vus_coordinates_category_10_set = set([val for sublist in brca_concurrent_vus_coordinates_category_10 for val in sublist])
    brca_concurrent_vus_coordinates_category_12_set = set([val for sublist in brca_concurrent_vus_coordinates_category_12 for val in sublist])
    brca_concurrent_vus_coordinates_category_13_set = set([val for sublist in brca_concurrent_vus_coordinates_category_13 for val in sublist])
    
    # Union all apparent benign vus variants into one set
    total_concurrent_vus_coordinates_set = set.union(brca_concurrent_vus_coordinates_category_8_set, brca_concurrent_vus_coordinates_category_10_set, brca_concurrent_vus_coordinates_category_13_set)
    
    # Union all samples
    
    ## Output co-occurrence report
    with open(options.outReport, 'w') as report_file:
        report_file.write("Pathogenic variant concordance\n")
        report_file.write("\t1|0 : {}\n".format(len(brca_pathogenic_left_het_sample_set)))
        report_file.write("\t0|1 : {}\n".format(len(brca_pathogenic_right_het_sample_set)))
        report_file.write("\t1|1 : {}\n".format(len(brca_pathogenic_hom_sample_set)))
        report_file.write("\n")
        report_file.write("VUS variant concordance\n")
        report_file.write("\t1|0 : {}\n".format(len(brca_vus_left_het_sample_set)))
        report_file.write("\t0|1 : {}\n".format(len(brca_vus_right_het_sample_set)))
        report_file.write("\t1|1 : {}\n".format(len(brca_vus_hom_sample_set)))
        report_file.write("\n")
        report_file.write("Pathogenic - VUS variant co-occurrence\n")
        report_file.write("\tPathogenic 1|0 - VUS 1|0 : {}\n".format(len(brca_left_het_coocourance_sample_set)))
        report_file.write("\tPathogenic 1|0 - VUS 0|1 : {}\n".format(len(brca_left_het_path_right_het_vus_coocourance_sample_set)))
        report_file.write("\tPathogenic 1|0 - VUS 1|1 : {}\n".format(len(brca_left_het_path_hom_vus_coocourance_sample_set)))
        report_file.write("\tPathogenic 0|1 - VUS 1|0 : {}\n".format(len(brca_right_het_path_left_het_vus_coocourance_sample_set)))
        report_file.write("\tPathogenic 0|1 - VUS 0|1 : {}\n".format(len(brca_right_het_path_right_het_vus_coocourance_sample_set)))
        report_file.write("\tPathogenic 0|1 - VUS 1|1 : {}\n".format(len(brca_right_het_path_hom_vus_coocourance_sample_set)))
        report_file.write("\tPathogenic 1|1 - VUS 1|0 : {}\n".format(len(brca_hom_path_left_het_vus_coocourance_sample_set)))
        report_file.write("\tPathogenic 1|1 - VUS 0|1 : {}\n".format(len(brca_hom_path_right_het_vus_coocourance_sample_set)))
        report_file.write("\tPathogenic 1|1 - VUS 1|1 : {}\n".format(len(brca_hom_path_hom_vus_coocourance_sample_set)))
        report_file.write("\n")
        report_file.write("Total unique samples in each trans category : {}\n".format(len(brca_trans_union_sample_set)))
        report_file.write("Total unique VUS in (Pathogenic 1|0 - VUS 0|1) set : {}\n".format(len(brca_concurrent_vus_coordinates_category_8_set)))
        report_file.write("Total unique VUS in (Pathogenic 1|0 - VUS 1|1) set : {}\n".format(len(brca_concurrent_vus_coordinates_category_9_set)))
        report_file.write("Total unique VUS in (Pathogenic 0|1 - VUS 1|0) set : {}\n".format(len(brca_concurrent_vus_coordinates_category_10_set)))
        report_file.write("Total unique VUS in (Pathogenic 0|1 - VUS 1|1) set : {}\n".format(len(brca_concurrent_vus_coordinates_category_12_set)))
        report_file.write("Total unique VUS in (VUS 1|1) set : {}\n".format(len(brca_concurrent_vus_coordinates_category_13_set)))
        report_file.write("Total unique apparent-benign VUS : {}\n".format(len(total_concurrent_vus_coordinates_set)))
    
    ## Output to verbose report:
    complete_report_filename = "complete_{}".format(options.outReport)
    with open(complete_report_filename, 'w') as complete_report_file:
        complete_report_file.write("sample_name\ttype\tapparent_benign_vus_variant\tsupporting_path_variants\n")
        # VUS variants in samples with (PATH 1|0 and VUS 0|1)
        for sample in brca_left_het_path_right_het_vus_coocourance_sample_set:
            complete_report_file.write("{}\tVUS_0|1_PATH_1|0\t{}\t{}\n".format(sample, brca_vus_right_het_sample_list[sample], brca_pathogenic_left_het_sample_list[sample]))
        # VUS variants in samples with (PATH 0|1 and VUS 1|0)
        for sample in brca_right_het_path_left_het_vus_coocourance_sample_set:
            complete_report_file.write("{}\tVUS_1|0_PATH_0|1\t{}\t{}\n".format(sample, brca_vus_left_het_sample_list[sample], brca_pathogenic_right_het_sample_list[sample]))
        # VUS variants in samples with (VUS 1|1)
        for sample in brca_vus_hom_sample_set:
            complete_report_file.write("{}\tVUS_1|1\t{}\n".format(sample, brca_vus_hom_sample_list[sample]))
    
        
    ## Output apparent-benign VUS variants list
    with open(options.outVariants, 'w') as vcf_file:
        vcf_file.write("##fileformat=VCFv4.2\n")
        if 'chr13' in record.CHROM:
            vcf_file.write("##contig=<ID=chr13,length=114364328>\n")
        elif 'chr17' in record.CHROM:
            vcf_file.write("##contig=<ID=chr17,length=83257441>\n")
        vcf_file.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate Allele Frequency from Best-guess Genotypes">\n')
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for variant in total_concurrent_vus_coordinates_set:
            variant_elements = variant.split('_')
            chromosome_element = variant_elements[0]
            position_element = variant_elements[1]
            ref_allele_element = variant_elements[2]
            alt_allele_element = variant_elements[3]
            allele_freq_element = variant_elements[4]
            vcf_file.write("{}\t{}\t.\t{}\t{}\t.\t.\t{}\n".format(chromosome_element,position_element,ref_allele_element,alt_allele_element,allele_freq_element))
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))

