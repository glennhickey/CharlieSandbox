import vcf, argparse, sys
from collections import defaultdict
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
    
    all_variant_hom_genotype_count = defaultdict(tuple)
    vus_variant_hom_genotype_count = defaultdict(tuple)
    vus_variant_lt_0_01_hom_genotype_count = defaultdict(tuple)
    vus_variant_lt_0_005_hom_genotype_count = defaultdict(tuple)
    vus_variant_lt_0_001_hom_genotype_count = defaultdict(tuple)
    num_samples = 0
    for record in vcf_reader:
        num_hom_alts = len(record.get_hom_alts())
        num_hom_refs = len(record.get_hom_refs())
        num_hets = len(record.get_hets())
        total_genotypes = float(num_hom_alts + num_hom_refs + num_hets)
        num_samples = total_genotypes
        hom_ref_frequency = float(num_hom_refs)/float(total_genotypes)
        hom_alt_frequency = float(num_hom_alts)/float(total_genotypes)
        het_frequency = float(num_hets)/float(total_genotypes)
        minor_allele_count = (num_hom_alts*2) + num_hets
        minor_af = float(minor_allele_count)/float((total_genotypes*2))
        variant_record = "{}_{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT[0],record.INFO['AF'][0])
        all_variant_hom_genotype_count[variant_record] = (num_hom_alts, minor_af)
        if variant_record in vus_list:
            vus_variant_hom_genotype_count[variant_record] = (num_hom_alts, minor_af)
            if minor_af <= 0.01:
                vus_variant_lt_0_01_hom_genotype_count[variant_record] = (num_hom_alts, minor_af)
            if minor_af <= 0.005:
                vus_variant_lt_0_005_hom_genotype_count[variant_record] = (num_hom_alts, minor_af)
            if minor_af <= 0.001:
                vus_variant_lt_0_001_hom_genotype_count[variant_record] = (num_hom_alts, minor_af)
    
    expected_hom_genotypes_all = 0.0
    expected_proportion_all = 0.0
    observed_hom_genotypes_all = 0.0
    for (hom_alt_count,minor_af) in all_variant_hom_genotype_count.values():
        observed_hom_genotypes_all += hom_alt_count
        expected_proportion_all += minor_af**2
    expected_hom_genotypes_all = expected_proportion_all*num_samples
    
    expected_hom_genotypes_vus = 0.0
    expected_proportion_vus = 0.0
    observed_hom_genotypes_vus = 0.0
    for (hom_alt_count,minor_af) in vus_variant_hom_genotype_count.values():
        observed_hom_genotypes_vus += hom_alt_count
        expected_proportion_vus += minor_af**2
    expected_hom_genotypes_vus = expected_proportion_vus*num_samples
    
    expected_hom_genotypes_vus_lt_0_01 = 0.0
    expected_proportion_vus_lt_0_01 = 0.0
    observed_hom_genotypes_vus_lt_0_01 = 0.0
    for (hom_alt_count,minor_af) in vus_variant_lt_0_01_hom_genotype_count.values():
        observed_hom_genotypes_vus_lt_0_01 += hom_alt_count
        expected_proportion_vus_lt_0_01 += minor_af**2
    expected_hom_genotypes_vus_lt_0_01 = expected_proportion_vus_lt_0_01*num_samples
    
    expected_hom_genotypes_vus_lt_0_005 = 0.0
    expected_proportion_vus_lt_0_005 = 0.0
    observed_hom_genotypes_vus_lt_0_005 = 0.0
    for (hom_alt_count,minor_af) in vus_variant_lt_0_005_hom_genotype_count.values():
        observed_hom_genotypes_vus_lt_0_005 += hom_alt_count
        expected_proportion_vus_lt_0_005 += minor_af**2
    expected_hom_genotypes_vus_lt_0_005 = expected_proportion_vus_lt_0_005*num_samples
    
    expected_hom_genotypes_vus_lt_0_001 = 0.0
    expected_proportion_vus_lt_0_001 = 0.0
    observed_hom_genotypes_vus_lt_0_001 = 0.0
    for (hom_alt_count,minor_af) in vus_variant_lt_0_001_hom_genotype_count.values():
        observed_hom_genotypes_vus_lt_0_001 += hom_alt_count
        expected_proportion_vus_lt_0_001 += minor_af**2
    expected_hom_genotypes_vus_lt_0_001 = expected_proportion_vus_lt_0_001*num_samples
    
    print('** ALL GENOTYPES **')
    print('expected number of homozygous alt genotypes: {}'.format(expected_hom_genotypes_all))
    print('observed number of homozygous alt genotypes: {}'.format(observed_hom_genotypes_all))
    
    print('** ALL VUS GENOTYPES **')
    print('expected number of homozygous alt genotypes: {}'.format(expected_hom_genotypes_vus))
    print('observed number of homozygous alt genotypes: {}'.format(observed_hom_genotypes_vus))
    
    print('** VUS WITH AF <= 0.01 GENOTYPES **')
    print('expected number of homozygous alt genotypes: {}'.format(expected_hom_genotypes_vus_lt_0_01))
    print('observed number of homozygous alt genotypes: {}'.format(observed_hom_genotypes_vus_lt_0_01))
    
    print('** VUS WITH AF <= 0.005 GENOTYPES **')
    print('expected number of homozygous alt genotypes: {}'.format(expected_hom_genotypes_vus_lt_0_005))
    print('observed number of homozygous alt genotypes: {}'.format(observed_hom_genotypes_vus_lt_0_005))
    
    print('** VUS WITH AF <= 0.001 GENOTYPES **')
    print('expected number of homozygous alt genotypes: {}'.format(expected_hom_genotypes_vus_lt_0_001))
    print('observed number of homozygous alt genotypes: {}'.format(observed_hom_genotypes_vus_lt_0_001))
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))

