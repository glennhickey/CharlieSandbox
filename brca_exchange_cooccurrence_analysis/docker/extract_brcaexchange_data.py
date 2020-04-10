from __future__ import print_function
import argparse, sys
from pprint import pprint  # for pretty-printing results
import requests as rq  # for issuing HTTP(S) queries

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input name of output brca exchange vcfs.')
    parser.add_argument('-o', '--outBrcaVCFnames', type=str,
        help='Output brcaexchange vcf base filename.')
    parser.add_argument('-g', '--gene', type=str,
        help='Extract which gene to process.')

    options = parser.parse_args()
    return options

def main(args):

    options = parse_args()
    
    contig_header_string="" 
    if options.gene == "BRCA1" :
        contig_header_string="##contig=<ID=chr17,length=83257441>\n"
    elif options.gene == "BRCA2" : 
        contig_header_string="##contig=<ID=chr13,length=114364328>\n"
        
    # Extract the brca pathogenic brcaexchange data and build the VCF
    brca_pathogenic_vcf_filename = "{}.{}.pathogenic.vcf".format(options.outBrcaVCFnames,options.gene)
    brca_pathogenic_query = """
    https://brcaexchange.org/backend/data/
    ?format=json
    &filter=Pathogenicity_expert
    &filterValue=Pathogenic
    &filter=Gene_Symbol
    &filterValue={}
    &order_by=Genomic_Coordinate_hg38
    &direction=ascending
    &search_term=
    &include=Variant_in_ENIGMA
    &include=Variant_in_ClinVar
    &include=Variant_in_1000_Genomes
    &include=Variant_in_ExAC&include=Variant_in_LOVD
    &include=Variant_in_BIC&include=Variant_in_ESP
    &include=Variant_in_exLOVD&include=Variant_in_Findlay_BRCA1_Ring_Function_Scores
    &include=Variant_in_GnomAD
    """.format(options.gene)
    with open(brca_pathogenic_vcf_filename, 'w') as vcf_file:
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write(contig_header_string)
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        brca_raw_data = rq.get(brca_pathogenic_query.replace(' ','').replace('\n', '')).json()
        for x in brca_raw_data['data']:
            vcf_file.write("chr{}\t{}\t.\t{}\t{}\t.\t.\t.\n".format(x['Chr'],x['Pos'],x['Ref'],x['Alt']))
    
    # Extract the brca variants of unknown significance brcaexchange data and build the VCF
    brca_vus_vcf_filename = "{}.{}.vus.vcf".format(options.outBrcaVCFnames,options.gene)
    brca_vus_query = """
    https://brcaexchange.org/backend/data/
    ?format=json
    &filter=Pathogenicity_expert
    &filterValue=Not_Yet_Reviewed
    &filter=Gene_Symbol
    &filterValue={}
    &order_by=Genomic_Coordinate_hg38
    &direction=ascending
    &search_term=
    &include=Variant_in_ENIGMA
    &include=Variant_in_ClinVar
    &include=Variant_in_1000_Genomes
    &include=Variant_in_ExAC&include=Variant_in_LOVD
    &include=Variant_in_BIC&include=Variant_in_ESP
    &include=Variant_in_exLOVD&include=Variant_in_Findlay_BRCA1_Ring_Function_Scores
    &include=Variant_in_GnomAD
    """.format(options.gene)
    with open(brca_vus_vcf_filename, 'w') as vcf_file:
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write(contig_header_string)
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        brca_raw_data = rq.get(brca_vus_query.replace(' ','').replace('\n', '')).json()
        for x in brca_raw_data['data']:
            vcf_file.write("chr{}\t{}\t.\t{}\t{}\t.\t.\t.\n".format(x['Chr'],x['Pos'],x['Ref'],x['Alt']))
    
    # Extract all of the brca variants brcaexchange data and build the VCF
    brca_all_vcf_filename = "{}.{}.all.vcf".format(options.outBrcaVCFnames,options.gene)
    brca_all_query = """
    https://brcaexchange.org/backend/data/
    ?format=json
    &filter=Gene_Symbol
    &filterValue={}
    &order_by=Genomic_Coordinate_hg38
    &direction=ascending
    &search_term=
    &include=Variant_in_ENIGMA
    &include=Variant_in_ClinVar
    &include=Variant_in_1000_Genomes
    &include=Variant_in_ExAC&include=Variant_in_LOVD
    &include=Variant_in_BIC&include=Variant_in_ESP
    &include=Variant_in_exLOVD&include=Variant_in_Findlay_BRCA1_Ring_Function_Scores
    &include=Variant_in_GnomAD
    """.format(options.gene)
    with open(brca_all_vcf_filename, 'w') as vcf_file:
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write(contig_header_string)
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        brca_raw_data = rq.get(brca_all_query.replace(' ','').replace('\n', '')).json()
        for x in brca_raw_data['data']:
            vcf_file.write("chr{}\t{}\t.\t{}\t{}\t.\t.\t.\n".format(x['Chr'],x['Pos'],x['Ref'],x['Alt']))
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))

