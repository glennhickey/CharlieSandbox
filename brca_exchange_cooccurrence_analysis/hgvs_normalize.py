import vcf, argparse, sys
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.location
import hgvs.edit, hgvs.posedit
import hgvs.sequencevariant
import hgvs.normalizer

import pyhgvs
import pyhgvs.utils as pyhgvs_utils
from pyfaidx import Fasta

from hgvs.extras.babelfish import Babelfish

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input a bgzip compressed and tabix indexed vcf and output hgvs normalized vcf filename.')
    parser.add_argument('-i', '--inVCF', type=str,
        help='Input bgzip compressed and tabix indexed vcf filepath.')
    parser.add_argument('-o', '--outVCF', type=str,
        help='Output hgvs-normalized VCF filename.')
    parser.add_argument('-r', '--refFASTA', type=str,
        help='Input FASTA format reference filename. ex: "hg38.p12.fa"')
    parser.add_argument('-g', '--refSEQ', type=str,
        help='Input GenePred format refSeq filename. ex: "ncbiRefSeq.txt"')

    options = parser.parse_args()
    return options

def main(args):

    options = parse_args()
    
    hdp = hgvs.dataproviders.uta.connect()
    am38 = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh38')
    hn = hgvs.normalizer.Normalizer(hdp)
    hp = hgvs.parser.Parser()
    # Read genome sequence using pyfaidx
    genome = Fasta(options.refFASTA)

    # Read RefSeq transcripts into a python dict.
    with open(options.refSEQ) as infile:
        transcripts = pyhgvs_utils.read_transcripts(infile)
    
    # Provide a callback for fetching a transcript by its name.
    def get_transcript(name):
        return transcripts.get(name)
    
    babelfish38 = Babelfish(hdp, assembly_name="GRCh38")
     
    ## extract base variant representation
    with open(options.inVCF, 'rb') as in_vcf, open(options.outVCF, 'w') as out_vcf:
        vcf_reader = vcf.Reader(in_vcf)
        vcf_writer = vcf.Writer(out_vcf, vcf_reader)
        for record in vcf_reader:
            # Convert variants for indel HGVS representation
            chrom, offset, ref, alt = (str(record.CHROM), record.POS, str(record.REF), str(record.ALT[0]))
            print('chrom: {}, offset: {}, ref: {}, alt: {}'.format(chrom, offset, ref, alt))
            if 'chr13' in record.CHROM:
                transcript_id = "NM_000059.3"
            elif 'chr17' in record.CHROM:
                transcript_id = "NM_007294.4"
            transcript = get_transcript(transcript_id)
            try:
                hgvs_name = pyhgvs.format_hgvs_name(chrom, offset, ref, alt, genome, transcript, use_gene=False, max_allele_length=50000)
                hgvs_c = hp.parse_hgvs_variant(hgvs_name)
                if len(ref) == len(alt) and len(ref) == 1:
                    # Variant is a SNP, normalize using hgvs Normalizer function
                    if 'chr17' in record.CHROM: hgvs_c.ac = 'NM_007294.3'
                    norm_hgvs_c = hn.normalize(hgvs_c)
                    if 'chr17' in record.CHROM: norm_hgvs_c.ac = 'NM_007294.4'
                    chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(str(norm_hgvs_c), genome, normalize=False, get_transcript=get_transcript)
                else:
                    # Variant is an INDEL, normalize using hgvs babelfish38.hgvs_to_vcf function
                    if 'chr17' in record.CHROM: hgvs_c.ac = 'NM_007294.3'
                    hgvs_g = am38.c_to_g(hgvs_c)
                    vcf_values = babelfish38.hgvs_to_vcf(hgvs_g)
                    chrom, offset, ref, alt = 'chr{}'.format(vcf_values[0]), vcf_values[1], vcf_values[2], vcf_values[3]
            except hgvs.exceptions.HGVSUnsupportedOperationError as e:
                print('hgvs.exceptions.HGVSUnsupportedOperationError: {}'.format(e))
            except hgvs.exceptions.HGVSInvalidIntervalError as e:
                print('hgvs.exceptions.HGVSInvalidIntervalError: {}'.format(e))
            except hgvs.exceptions.HGVSInvalidVariantError as e:
                print('hgvs.exceptions.HGVSInvalidVariantError: {}'.format(e))
            except AttributeError as e:
                print('AttributeError: {}'.format(e))
            except KeyError as e:
                print('KeyError: {}'.format(e))
            # Update and write the new normalized record
            record.POS = offset
            record.REF = ref
            record.ALT = [alt]
            vcf_writer.write_record(record)
            #print('')
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))

