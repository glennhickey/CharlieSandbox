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
##      - rtg is present in $PROGRAM_DIR
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
        rtg                             (from https://github.com/RealTimeGenomics/rtg-tools)
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

######################################
## STRATIFY VCFS BY SNPS AND INDELS ##
echo "stratifying vcfs by snps and indels"
split_vcf () {
    VCF_FILE="$1"
    OUTPUT_BASE_FILENAME="$2"
    PROGRAM_DIR="$3"
    ${PROGRAM_DIR}/bcftools view --types snps ${VCF_FILE} > ${OUTPUT_BASE_FILENAME}.snps.vcf
    bgzip ${OUTPUT_BASE_FILENAME}.snps.vcf
    tabix -p vcf ${OUTPUT_BASE_FILENAME}.snps.vcf.gz
    ${PROGRAM_DIR}/bcftools view --types indels ${VCF_FILE} > ${OUTPUT_BASE_FILENAME}.indels.vcf
    bgzip ${OUTPUT_BASE_FILENAME}.indels.vcf
    tabix -p vcf ${OUTPUT_BASE_FILENAME}.indels.vcf.gz
}

# stratify BRCA1 brcaexchange data
split_vcf brcaexchange_variants.brca1.pathogenic.vcf.gz "brcaexchange_variants.brca1.pathogenic" ${PROGRAM_DIR}
split_vcf brcaexchange_variants.brca1.vus.vcf.gz "brcaexchange_variants.brca1.vus" ${PROGRAM_DIR}

# stratify BRCA2 brcaexchange data
split_vcf brcaexchange_variants.brca2.pathogenic.vcf.gz "brcaexchange_variants.brca2.pathogenic" ${PROGRAM_DIR}
split_vcf brcaexchange_variants.brca2.vus.vcf.gz "brcaexchange_variants.brca2.vus" ${PROGRAM_DIR}

# stratify BRCA1 sample set data
split_vcf raw_sample.phased.brca1region.vcf.gz "raw_sample.phased.brca1region" ${PROGRAM_DIR}

# stratify BRCA12 sample set data
split_vcf raw_sample.phased.brca2region.vcf.gz "raw_sample.phased.brca2region" ${PROGRAM_DIR}
echo "stratifying vcfs by snps and indels DONE"
echo ""

################################
## RUN VCF COMPARISON FILTERS ##
echo "running vcf comparison filters"
run_rtg_vcfeval () {
    BASE_VCF="$1"
    QUERY_VCF="$2"
    REF_INDEX="$3"
    RTG_VCFEVAL_WORKDIR="$4"
    OUTPUT_FILENAME="$5"
    PROGRAM_DIR="$6"
    rmdir ${RTG_VCFEVAL_WORKDIR}
    ${PROGRAM_DIR}/rtg vcfeval \
        -T 16 \
        --output-mode=split \
        --sample ALT \
        --ref-overlap \
        --XXcom.rtg.vcf.eval.max-paths=500000 \
        --XXcom.rtg.vcf.eval.max-iterations=10000000 \
         -b ${BASE_VCF} \
         -c ${QUERY_VCF} \
         -t ${REF_INDEX} \
         -o ${RTG_VCFEVAL_WORKDIR}/vcfeval_output
    ${PROGRAM_DIR}/vcf-sort ${RTG_VCFEVAL_WORKDIR}/vcfeval_output/tp.vcf.gz > ${OUTPUT_FILENAME}
}

# Process brca1 pathogenic variants
BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca1.pathogenic.snps.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca1region.snps.vcf.gz"
REF_INDEX="${WORK_DIR}/hg38.p12.fa.sdf"
RTG_VCFEVAL_WORKDIR="${WORK_DIR}/brca1_PATHOGENIC_SNPS"
OUTPUT_FILENAME="brca1_PATHOGENIC_SNPS.sorted.vcf"
run_rtg_vcfeval ${BASE_VCF} ${QUERY_VCF} ${REF_INDEX} ${RTG_VCFEVAL_WORKDIR} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca1.pathogenic.indels.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca1region.indels.vcf.gz"
RTG_VCFEVAL_WORKDIR="${WORK_DIR}/brca1_PATHOGENIC_INDELS"
OUTPUT_FILENAME="brca1_PATHOGENIC_INDELS.sorted.vcf"
run_rtg_vcfeval ${BASE_VCF} ${QUERY_VCF} ${REF_INDEX} ${RTG_VCFEVAL_WORKDIR} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

# Process brca1 vus variants
BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca1.vus.snps.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca1region.snps.vcf.gz"
RTG_VCFEVAL_WORKDIR="${WORK_DIR}/brca1_VUS_SNPS"
OUTPUT_FILENAME="brca1_VUS_SNPS.sorted.vcf"
run_rtg_vcfeval ${BASE_VCF} ${QUERY_VCF} ${REF_INDEX} ${RTG_VCFEVAL_WORKDIR} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca1.vus.indels.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca1region.indels.vcf.gz"
RTG_VCFEVAL_WORKDIR="${WORK_DIR}/brca1_VUS_INDELS"
OUTPUT_FILENAME="brca1_VUS_INDELS.sorted.vcf"
run_rtg_vcfeval ${BASE_VCF} ${QUERY_VCF} ${REF_INDEX} ${RTG_VCFEVAL_WORKDIR} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

# Process brca2 pathogenic variants
BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca2.pathogenic.snps.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca2region.snps.vcf.gz"
RTG_VCFEVAL_WORKDIR="${WORK_DIR}/brca2_PATHOGENIC_SNPS"
OUTPUT_FILENAME="brca2_PATHOGENIC_SNPS.sorted.vcf"
run_rtg_vcfeval ${BASE_VCF} ${QUERY_VCF} ${REF_INDEX} ${RTG_VCFEVAL_WORKDIR} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca2.pathogenic.indels.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca2region.indels.vcf.gz"
RTG_VCFEVAL_WORKDIR="${WORK_DIR}/brca2_PATHOGENIC_INDELS"
OUTPUT_FILENAME="brca2_PATHOGENIC_INDELS.sorted.vcf"
run_rtg_vcfeval ${BASE_VCF} ${QUERY_VCF} ${REF_INDEX} ${RTG_VCFEVAL_WORKDIR} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

# Process brca2 vus variants
BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca2.vus.snps.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca2region.snps.vcf.gz"
RTG_VCFEVAL_WORKDIR="${WORK_DIR}/brca2_VUS_SNPS"
OUTPUT_FILENAME="brca2_VUS_SNPS.sorted.vcf"
run_rtg_vcfeval ${BASE_VCF} ${QUERY_VCF} ${REF_INDEX} ${RTG_VCFEVAL_WORKDIR} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

BASE_VCF="${WORK_DIR}/brcaexchange_variants.brca2.vus.indels.vcf.gz"
QUERY_VCF="${WORK_DIR}/raw_sample.phased.brca2region.indels.vcf.gz"
RTG_VCFEVAL_WORKDIR="${WORK_DIR}/brca2_VUS_INDELS"
OUTPUT_FILENAME="brca2_VUS_INDELS.sorted.vcf"
run_rtg_vcfeval ${BASE_VCF} ${QUERY_VCF} ${REF_INDEX} ${RTG_VCFEVAL_WORKDIR} ${OUTPUT_FILENAME} ${PROGRAM_DIR}
echo "running vcf comparison filters DONE"
echo ""

####################################
## CONCATENATE SNP AND INDEL VCFS ##
echo "concatenating snp and indel vcfs"
merge_vcfs () {
    SNPS_VCF="$1"
    INDELS_VCF="$2"
    OUTPUT_BASE_FILENAME="$3"
    PROGRAM_DIR="$4"
    echo "$SNPS_VCF"
    echo "$INDELS_VCF"
    echo "$OUTPUT_BASE_FILENAME"
    echo "$PROGRAM_DIR"
    ${PROGRAM_DIR}/bcftools concat ${SNPS_VCF} ${INDELS_VCF} > ${OUTPUT_BASE_FILENAME}.vcf
    ${PROGRAM_DIR}/vcf-sort ${OUTPUT_BASE_FILENAME}.vcf > ${OUTPUT_BASE_FILENAME}.sorted.vcf
    bgzip ${OUTPUT_BASE_FILENAME}.sorted.vcf
    tabix -p vcf ${OUTPUT_BASE_FILENAME}.sorted.vcf.gz
    rm -f ${OUTPUT_BASE_FILENAME}.vcf
    echo "${OUTPUT_BASE_FILENAME}.sorted.vcf.gz"
}

# merge and sort matching snps and indels from brca1 pathogenic data
brca1_PATHOGENIC_SNPS_file="brca1_PATHOGENIC_SNPS.sorted.vcf"
brca1_PATHOGENIC_INDELS_file="brca1_PATHOGENIC_INDELS.sorted.vcf"
OUTPUT_BASE_FILENAME="brca1_PATHOGENIC_SNPS_INDELS"
merge_vcfs ${brca1_PATHOGENIC_SNPS_file} ${brca1_PATHOGENIC_INDELS_file} ${OUTPUT_BASE_FILENAME} ${PROGRAM_DIR}

# merge and sort matching snps and indels from brca1 vus data
brca1_VUS_SNPS_file="brca1_VUS_SNPS.sorted.vcf"
brca1_VUS_INDELS_file="brca1_VUS_INDELS.sorted.vcf"
OUTPUT_BASE_FILENAME="brca1_VUS_SNPS_INDELS"
merge_vcfs ${brca1_VUS_SNPS_file} ${brca1_VUS_INDELS_file} ${OUTPUT_BASE_FILENAME} ${PROGRAM_DIR}

# merge and sort matching snps and indels from brca2 pathogenic data
brca2_PATHOGENIC_SNPS_file="brca2_PATHOGENIC_SNPS.sorted.vcf"
brca2_PATHOGENIC_INDELS_file="brca2_PATHOGENIC_INDELS.sorted.vcf"
OUTPUT_BASE_FILENAME="brca2_PATHOGENIC_SNPS_INDELS"
merge_vcfs ${brca2_PATHOGENIC_SNPS_file} ${brca2_PATHOGENIC_INDELS_file} ${OUTPUT_BASE_FILENAME} ${PROGRAM_DIR}

# merge and sort matching snps and indels from brca2 vus data
brca2_VUS_SNPS_file="brca2_VUS_SNPS.sorted.vcf"
brca2_VUS_INDELS_file="brca2_VUS_INDELS.sorted.vcf"
OUTPUT_BASE_FILENAME="brca2_VUS_SNPS_INDELS"
merge_vcfs ${brca2_VUS_SNPS_file} ${brca2_VUS_INDELS_file} ${OUTPUT_BASE_FILENAME} ${PROGRAM_DIR}
echo "concatenating snp and indel vcfs DONE"
echo ""

##################################
## RUN THE CONCORDANCE ANALYSIS ##
echo "running the concordance analysis"
source ${WORK_DIR}/vcfvenv/bin/activate
python ${PROGRAM_DIR}/detect_vus_benign.py -i "${WORK_DIR}/brca1_VUS_SNPS_INDELS.sorted.vcf.gz" -j "${WORK_DIR}/brca1_PATHOGENIC_SNPS_INDELS.sorted.vcf.gz" -o "${WORK_DIR}/brca1_coocurrence_report.txt"
python ${PROGRAM_DIR}/detect_vus_benign.py -i "${WORK_DIR}/brca2_VUS_SNPS_INDELS.sorted.vcf.gz" -j "${WORK_DIR}/brca2_PATHOGENIC_SNPS_INDELS.sorted.vcf.gz" -o "${WORK_DIR}/brca2_coocurrence_report.txt"
deactivate
echo "running the concordance analysis DONE"
echo ""

exit

