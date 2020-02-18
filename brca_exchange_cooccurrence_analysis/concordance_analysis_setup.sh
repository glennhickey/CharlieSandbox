#!/bin/bash
#################################################################################################
##
##  Script to setup and run vcf filtering on phased population genotype VCF data prior to
##	running VUS - PATHOGENIC variant concordance analysis.
##
##  Inputs:
##
##  Assumptions:
##      - Everything is in GRCh38-based coordinates
##      - Raw genotype data is phased
##      - bgzip is installed and accessible in the $PATH environment variable
##      - tabix is installed and accessible in the $PATH environment variable
##      - extract_brcaexchange_data.py is present in $PROGRAM_DIR
##      - detect_vus_benign.py is present in $PROGRAM_DIR
##      - bcftools is present in $PROGRAM_DIR
##      - vcf-sort is present in $PROGRAM_DIR
##
##  Last modified:
##  Last modified by: Charles Markello
##
#################################################################################################

## Create help statement
usage(){
cat << EOF
This script setups up input directories, downloads files, and runs VUS - Pathogenic concordance analysis.
Inputs:
    -a PATH to raw population joint-genotyped chromosome 17 VCF file in BCF format.
    -b PATH to raw population joint-genotyped chromosome 13 VCF file in BCF format.
    -p PATH to where programs and scripts live. Must have the following:
        extract_brcaexchange_data.py    (from https://github.com/cmarkello/Sandbox/tree/master/brca_exchange_cooccurrence_analysis)
        detect_vus_benign.py            (from https://github.com/cmarkello/Sandbox/tree/master/brca_exchange_cooccurrence_analysis)
        bcftools                        (from http://www.htslib.org/download/)
        vcf-sort                        (from https://vcftools.github.io/perl_module.html)
    -w PATH to where the working directory will be where data will be processed.
    
Outputs:
Assumptions:
EOF

}

## Check number of arguments
if [ $# -lt 4 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## Parse through arguments
while getopts "a:b:p:w:h" OPTION; do
    case $OPTION in
        a)
            CHR17_BCF=$OPTARG
        ;;
        b)
            CHR13_BCF=$OPTARG
        ;;
        p)
            PROGRAM_DIR=$OPTARG
        ;;
        w)
            WORK_DIR=$OPTARG
        ;;
        h)
            usage
            exit 1
        ;;
        \?)
            usage
            exit 1
        ;;
    esac
done

if [ ! -d "${PROGRAM_DIR}" ]; then
    mkdir -p ${PROGRAM_DIR}
    chmod 2770 ${PROGRAM_DIR}
fi

if [ ! -d "${WORK_DIR}" ]; then
    mkdir -p ${WORK_DIR}
    chmod 2770 ${WORK_DIR}
fi
cd ${WORK_DIR}

########################################################################
## build the python virtualenvironment and install necessary packages ##
echo "building python virtualenv"
virtualenv ${WORK_DIR}/vcfvenv
source ${WORK_DIR}/vcfvenv/bin/activate
pip install requests pyvcf
deactivate
echo "building python virtualenv DONE"
echo ""

##############################################
## download brcaexchange data to vcf format ##
echo "downloading and formatting brcaexchange data"
source ${WORK_DIR}/vcfvenv/bin/activate
python ${PROGRAM_DIR}/extract_brcaexchange_data.py -o brcaexchange_variants
deactivate
# Create the brcaexchange brca1 pathogenic and VUS list vcf files
bgzip brcaexchange_variants.brca1.pathogenic.vcf
tabix -p vcf brcaexchange_variants.brca1.pathogenic.vcf.gz
bgzip brcaexchange_variants.brca1.vus.vcf
tabix -p vcf brcaexchange_variants.brca1.vus.vcf.gz

# Create the brcaexchange brca2 pathogenic and VUS list vcf files
bgzip brcaexchange_variants.brca2.pathogenic.vcf
tabix -p vcf brcaexchange_variants.brca2.pathogenic.vcf.gz
bgzip brcaexchange_variants.brca2.vus.vcf
tabix -p vcf brcaexchange_variants.brca2.vus.vcf.gz
echo "downloading and formatting brcaexchange data DONE"
echo ""

############################
## FORMAT REFERENCE INDEX ##
echo "formatting reference index"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.fa.gz .
gzip -df hg38.p12.fa.gz > hg38.p12.fa
${PROGRAM_DIR}/rtg format -o hg38.p12.fa.sdf hg38.p12.fa
echo "formatting reference index DONE"
echo ""

#####################################
## FORMAT RAW SAMPLE GENOTYPE VCFS ##
echo "formatting raw sample genotype vcfs"
# Extract brca1 region (GRCh38 coordinates)
${PROGRAM_DIR}/bcftools filter --regions chr17:43044295-43170245 ${CHR17_BCF} > raw_sample.phased.brca1region.vcf
bgzip raw_sample.phased.brca1region.vcf
tabix -p vcf raw_sample.phased.brca1region.vcf.gz

# Extract brca2 region (GRCh38 coordinates)
${PROGRAM_DIR}/bcftools filter --regions chr13:32310000-32400000 ${CHR13_BCF} > raw_sample.phased.brca2region.vcf
bgzip raw_sample.phased.brca2region.vcf
tabix -p vcf raw_sample.phased.brca2region.vcf.gz
echo "formatting raw sample genotype vcfs DONE"
echo ""

################################
## RUN VCF COMPARISON FILTERS ##
echo "running vcf comparison filters"
run_vcf_comparison () {
    BASE_VCF="$1"
    QUERY_VCF="$2"
    OUTPUT_FILENAME="$3"
    PROGRAM_DIR="$4"
    ${PROGRAM_DIR}/bcftools isec \
        -O v \
        -n =2 -w 1 \
        ${QUERY_VCF} \
        ${BASE_VCF} \
        > ${OUTPUT_FILENAME}.vcf
    ${PROGRAM_DIR}/vcf-sort ${OUTPUT_FILENAME}.vcf > ${OUTPUT_FILENAME}.sorted.vcf
    bgzip ${OUTPUT_FILENAME}.sorted.vcf
    tabix -p vcf ${OUTPUT_FILENAME}.sorted.vcf.gz
    rm ${OUTPUT_FILENAME}.vcf
}

# Process brca1 pathogenic variants
BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca1.pathogenic.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca1region.vcf.gz"
OUTPUT_FILENAME="brca1_PATHOGENIC"
run_vcf_comparison ${BASE_VCF} ${QUERY_VCF} ${OUTPUT_FILENAME} ${PROGRAM_DIR}


# Process brca1 vus variants
BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca1.vus.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca1region.vcf.gz"
OUTPUT_FILENAME="brca1_VUS"
run_vcf_comparison ${BASE_VCF} ${QUERY_VCF} ${OUTPUT_FILENAME} ${PROGRAM_DIR}


# Process brca2 pathogenic variants
BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca2.pathogenic.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca2region.vcf.gz"
OUTPUT_FILENAME="brca2_PATHOGENIC"
run_vcf_comparison ${BASE_VCF} ${QUERY_VCF} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

# Process brca2 vus variants
BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca2.vus.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca2region.vcf.gz"
OUTPUT_FILENAME="brca2_VUS"
run_vcf_comparison ${BASE_VCF} ${QUERY_VCF} ${OUTPUT_FILENAME} ${PROGRAM_DIR}
echo "running vcf comparison filters DONE"
echo ""

##################################
## RUN THE CONCORDANCE ANALYSIS ##
echo "running the concordance analysis"
source ${WORK_DIR}/vcfvenv/bin/activate
python ${PROGRAM_DIR}/detect_vus_benign.py -i "${WORK_DIR}/brca1_VUS.sorted.vcf.gz" -j "${WORK_DIR}/brca1_PATHOGENIC.sorted.vcf.gz" -o "${WORK_DIR}/brca1_coocurrence_report.txt" -v "${WORK_DIR}/brca1_apparent_benign_VUS.vcf"
python ${PROGRAM_DIR}/detect_vus_benign.py -i "${WORK_DIR}/brca2_VUS.sorted.vcf.gz" -j "${WORK_DIR}/brca2_PATHOGENIC.sorted.vcf.gz" -o "${WORK_DIR}/brca2_coocurrence_report.txt" -v "${WORK_DIR}/brca2_apparent_benign_VUS.vcf"
deactivate
${PROGRAM_DIR}/vcf-sort ${WORK_DIR}/brca1_apparent_benign_VUS.vcf > ${WORK_DIR}/brca1_apparent_benign_VUS.sorted.vcf
bgzip ${WORK_DIR}/brca1_apparent_benign_VUS.sorted.vcf
tabix -p vcf ${WORK_DIR}/brca1_apparent_benign_VUS.sorted.vcf.gz
rm ${WORK_DIR}/brca1_apparent_benign_VUS.vcf
${PROGRAM_DIR}/vcf-sort ${WORK_DIR}/brca2_apparent_benign_VUS.vcf > ${WORK_DIR}/brca2_apparent_benign_VUS.sorted.vcf
bgzip ${WORK_DIR}/brca2_apparent_benign_VUS.sorted.vcf
tabix -p vcf ${WORK_DIR}/brca2_apparent_benign_VUS.sorted.vcf.gz
rm ${WORK_DIR}/brca2_apparent_benign_VUS.vcf
echo "running the concordance analysis DONE"
echo ""

exit

