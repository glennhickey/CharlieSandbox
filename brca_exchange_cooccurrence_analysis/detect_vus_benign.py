import vcf, argparse, sys
from collections import defaultdict

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input both a VUS and pathogenic concordant bgzip compressed and tabix indexed vcf.')
    parser.add_argument('-i', '--inVUSvcf', type=argparse.FileType('r'),
        help='Input brcaexchange-VUS/sample-genotype concordant vcf filepath.')
    parser.add_argument('-j', '--inPATHvcf', type=argparse.FileType('r'),
        help='Input brcaexchange-pathogenic/sample-genotype concordant vcf filepath.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')
    parser.add_argument('-v', '--outVariants', type=str,
        help='Output apparent-benign VUS variants list filename.')

    options = parser.parse_args()
    return options

def main(args):

    options = parse_args()

    ## Isolate samples by the 2 different phased PATHOGENIC het calls
    vcf_reader_pathogenic = vcf.Reader(options.inPATHvcf)
    brca_pathogenic_left_het_sample_list = defaultdict(list)
    brca_pathogenic_right_het_sample_list = defaultdict(list)
    brca_pathogenic_hom_sample_list = defaultdict(list)
    for record in vcf_reader_pathogenic:
        for sample in record.samples:
            if sample['GT'] == '1|0':
                brca_pathogenic_left_het_sample_list[sample.sample].append("{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT))
            elif sample['GT'] == '0|1':
                brca_pathogenic_right_het_sample_list[sample.sample].append("{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT))
            elif sample['GT'] == '1|1':
                brca_pathogenic_hom_sample_list[sample.sample].append("{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT))
    
    ## Isolate samples by the 2 different phased VUS het calls
    vcf_reader_vus = vcf.Reader(options.inVUSvcf)
    brca_vus_left_het_sample_list = defaultdict(list)
    brca_vus_right_het_sample_list = defaultdict(list)
    brca_vus_hom_sample_list = defaultdict(list)
    for record in vcf_reader_vus:
        for sample in record.samples:
            if sample['GT'] == '1|0':
                brca_vus_left_het_sample_list[sample.sample].append("{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT))
            elif sample['GT'] == '0|1':
                brca_vus_right_het_sample_list[sample.sample].append("{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT))
            elif sample['GT'] == '1|1':
                brca_vus_hom_sample_list[sample.sample].append("{}_{}_{}_{}".format(record.CHROM,record.POS,record.REF,record.ALT))

    ## Look for shared samples between cis het VUS and PATHOGENIC variants
    brca_pathogenic_left_het_sample_set = set(brca_pathogenic_left_het_sample_list.keys())
    brca_pathogenic_right_het_sample_set = set(brca_pathogenic_right_het_sample_list.keys())
    brca_pathogenic_hom_sample_set = set(brca_pathogenic_hom_sample_list.keys())
    brca_vus_left_het_sample_set = set(brca_vus_left_het_sample_list.keys())
    brca_vus_right_het_sample_set = set(brca_vus_right_het_sample_list.keys())
    brca_vus_hom_sample_set = set(brca_vus_hom_sample_list.keys())

    brca_left_het_path_right_het_path_set = brca_pathogenic_left_het_sample_set.intersection(brca_pathogenic_right_het_sample_set)
    brca_left_het_coocourance_sample_set = brca_pathogenic_left_het_sample_set.intersection(brca_vus_left_het_sample_set)
    brca_left_het_path_right_het_vus_coocourance_sample_set = brca_pathogenic_left_het_sample_set.intersection(brca_vus_right_het_sample_set)
    brca_left_het_path_hom_vus_coocourance_sample_set = brca_pathogenic_left_het_sample_set.intersection(brca_vus_hom_sample_set)
    brca_right_het_path_left_het_vus_coocourance_sample_set = brca_pathogenic_right_het_sample_set.intersection(brca_vus_left_het_sample_set)
    brca_right_het_path_right_het_vus_coocourance_sample_set = brca_pathogenic_right_het_sample_set.intersection(brca_vus_right_het_sample_set)
    brca_right_het_path_hom_vus_coocourance_sample_set = brca_pathogenic_right_het_sample_set.intersection(brca_vus_hom_sample_set)
    brca_hom_path_left_het_vus_coocourance_sample_set = brca_pathogenic_hom_sample_set.intersection(brca_vus_left_het_sample_set)
    brca_hom_path_right_het_vus_coocourance_sample_set = brca_pathogenic_hom_sample_set.intersection(brca_vus_right_het_sample_set)
    brca_hom_path_hom_vus_coocourance_sample_set = brca_pathogenic_hom_sample_set.intersection(brca_vus_hom_sample_set)

    brca_trans_intersection_sample_set = set.intersection(brca_left_het_path_right_het_vus_coocourance_sample_set, brca_right_het_path_left_het_vus_coocourance_sample_set, brca_right_het_path_hom_vus_coocourance_sample_set, brca_left_het_path_hom_vus_coocourance_sample_set)

    brca_trans_union_sample_set = set.union(brca_left_het_path_right_het_vus_coocourance_sample_set, brca_right_het_path_left_het_vus_coocourance_sample_set, brca_right_het_path_hom_vus_coocourance_sample_set, brca_left_het_path_hom_vus_coocourance_sample_set)

    ## Look at unique vus variants that coincide with concurrance sample set
    brca_concurrent_vus_coordinates_category_8 = list()
    for sample in brca_left_het_path_right_het_vus_coocourance_sample_set:
        brca_concurrent_vus_coordinates_category_8.append(brca_vus_right_het_sample_list[sample])

    brca_concurrent_vus_coordinates_category_9 = list()
    for sample in brca_left_het_path_hom_vus_coocourance_sample_set:
        brca_concurrent_vus_coordinates_category_9.append(brca_vus_hom_sample_list[sample])

    brca_concurrent_vus_coordinates_category_10 = list()
    for sample in brca_right_het_path_left_het_vus_coocourance_sample_set:
        brca_concurrent_vus_coordinates_category_10.append(brca_vus_left_het_sample_list[sample])


    brca_concurrent_vus_coordinates_category_12 = list()
    for sample in brca_right_het_path_hom_vus_coocourance_sample_set:
        brca_concurrent_vus_coordinates_category_12.append(brca_vus_hom_sample_list[sample])


    brca_concurrent_vus_coordinates_category_8_set = set([val for sublist in brca_concurrent_vus_coordinates_category_8 for val in sublist])
    brca_concurrent_vus_coordinates_category_9_set = set([val for sublist in brca_concurrent_vus_coordinates_category_9 for val in sublist])
    brca_concurrent_vus_coordinates_category_10_set = set([val for sublist in brca_concurrent_vus_coordinates_category_10 for val in sublist])
    brca_concurrent_vus_coordinates_category_12_set = set([val for sublist in brca_concurrent_vus_coordinates_category_12 for val in sublist])

    total_concurrent_vus_coordinates_set = set.union(brca_concurrent_vus_coordinates_category_8_set, brca_concurrent_vus_coordinates_category_9_set, brca_concurrent_vus_coordinates_category_10_set, brca_concurrent_vus_coordinates_category_12_set)

    # total_concurrent_vus_coordinates_set represents the vus brca exchange variants where at least one sample in the topMed dataset had that vus variant called a het at one spot and also had a pathogenic brca exchange variant that was trans het at another spot.

    ## Count number of unique apparent benign VUS variants that coocour with > 1 pathogenic variant
    #brca_aparent_benign_VUS_coocour_greater_1_PATH_set = set()
    #for sample in brca_left_het_path_right_het_vus_coocourance_sample_set:
    #    if len(brca_pathogenic_left_het_sample_list[sample]) > 1:
    #        brca_aparent_benign_VUS_coocour_greater_1_PATH_set.append()

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
        report_file.write("Total unique apparent-benign VUS : {}\n".format(len(total_concurrent_vus_coordinates_set)))
        
    ## Output apparent-benign VUS variants list
    with open(options.outVariants, 'w') as vcf_file:
        vcf_file.write("##fileformat=VCFv4.2\n")
        if 'chr13' in record.CHROM:
            vcf_file.write("##contig=<ID=chr13,length=114364328>\n")
        elif 'chr17' in record.CHROM:
            vcf_file.write("##contig=<ID=chr17,length=83257441>\n")
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for variant in total_concurrent_vus_coordinates_set:
            variant_elements = variant.split('_')
            chromosome_element = variant_elements[0]
            position_element = variant_elements[1]
            ref_allele_element = variant_elements[2]
            alt_allele_element = variant_elements[3].strip('][').split(', ')[0]
            vcf_file.write("{}\t{}\t.\t{}\t{}\t.\t.\t.\n".format(chromosome_element,position_element,ref_allele_element,alt_allele_element))
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))

