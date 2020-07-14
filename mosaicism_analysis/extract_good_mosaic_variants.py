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
    parser.add_argument('-s', '--sibling_names', nargs='?', type=str, default=[], action='append',
        help="sample names of siblings, including the proband. Optional.")
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')

    options = parser.parse_args()
    return options

def main(args):
    options = parse_args()

    vcf_reader = vcf.Reader(open(options.inVCF, 'rb'))
    # Collect sibling sample names
    sibling_name_list = list()
    if len(options.sibling_names) >= 1:
        sibling_name_list = options.sibling_names
    else:
        for sample_name in vcf_reader.samples:
            if sample_name not in [options.inMaternalID, options.inPaternalID]:
                sibling_name_list.append(sample_name)
    
    sibling_records = defaultdict(list)
    sibling_filestreams = {}
    for sibling_name in sibling_name_list:
        sibling_filestreams[sibling_name] = open("{}.{}.tsv".format(sibling_name,options.outReport), 'wb')
        sibling_filestreams[sibling_name].write(b'#CHR\tPOS\tMATERNAL_AD\tPATERNAL_AD\n')
        
    for record in vcf_reader:
        # Only include sites if it is a snp with only one alternative allele
        if not record.is_snp or not len(record.ALT) == 1: continue
        maternal_call = record.genotype(options.inMaternalID)
        paternal_call = record.genotype(options.inPaternalID)
        # Only include records for any siblings if both parents have genotype calls
        if not maternal_call.called or not paternal_call.called: continue
        # Only include records for any siblings if both parents total read depths are between 20 and 150
        if (min([sum(paternal_call.data.AD), sum(maternal_call.data.AD)]) < 20) or (max([sum(paternal_call.data.AD), sum(maternal_call.data.AD)]) > 150): continue
        # Only include records for any siblings if both parents have phasable genotypes
        #   Either both are opposite homozygous (0/0 and 1/1) or (1/1 and 0/0)
        #   If only one parent is heterozygous and the other parent is homozygous ([0/0 or 1/1] and [0/1 or 1/0]) or ([0/1 or 1/0] and [0/0 or 1/1])
        if not ((maternal_call.is_het and not paternal_call.is_het) or (paternal_call.is_het and not maternal_call.is_het)): continue
        for sibling_name in sibling_name_list:
            sibling_call = record.genotype(sibling_name)
            # Only include records for a sibling if it actually has a genotype call
            if not sibling_call.called: continue
            # Only include records for a sibling if the genotype called is a het
            if not sibling_call.is_het: continue
            # Only include records for a sibling if its total read depth is between 20 and 150
            if (sum(sibling_call.data.AD) < 20) or (sum(sibling_call.data.AD) > 150): continue
            # Extract phased read depth of the sibling
            if (not maternal_call.is_het and sibling_call.data.GT[0] == maternal_call.data.GT[0]) or (not paternal_call.is_het and sibling_call.data.GT[2] == paternal_call.data.GT[2]):
                sample_maternal_phased_AD = sibling_call.data.AD[0]
                sample_paternal_phased_AD = sibling_call.data.AD[1]
            else:
                sample_maternal_phased_AD = sibling_call.data.AD[1]
                sample_paternal_phased_AD = sibling_call.data.AD[0]
            sibling_filestreams[sibling_name].write("{}\t{}\t{}\t{}\n".format(record.CHROM, record.POS, sample_maternal_phased_AD, sample_paternal_phased_AD).encode())
    
    for sibling_name in sibling_name_list:
        sibling_filestreams[sibling_name].close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

