import vcf, argparse, sys
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
    parser.add_argument('-m', '--inMaternalID', type=str,
        help='Input maternal vcf sample name.')
    parser.add_argument('-p', '--inPaternalID', type=str,
        help='Inputpmaternal vcf sample name.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')

    options = parser.parse_args()
    return options

def main(args):
    
    options = parse_args()

    vcf_reader = vcf.Reader(open(options.inVCF, 'rb'))
    # Collect sibling sample names
    sibling_name_list = list()
    for sample_name in vcf_reader.samples:
        if sample_name not in [options.inMaternalID, options.inPaternalID]
            sibling_name_list.append(sample_name)
    
    with open(options.outReport, 'w') as report_file:
        for record in vcf_reader:
            if not (record.is_snp or len(record.ALT) == 1): continue
            maternal_call = record.genotype(options.inMaternalID)
            paternal_call = record.genotype(options.inPaternalID)
            for sibling_name in sibling_name_list:
                sibling_call = record.genotype(sibling_name)
                
        report_file.write("Number of samples: {}\n".format(num_samples))
        report_file.write("Number of homozygous alt genotypes for a sample:\n")
        report_file.write("Min: {}\n".format(min_genotypes))
        report_file.write("1st Quartile: {}\n".format(Q1_genotypes))
        report_file.write("Median: {}\n".format(median_genotypes))
        report_file.write("3rd Quartile: {}\n".format(Q3_genotypes))
        report_file.write("Max: {}\n".format(max_genotypes))
    

