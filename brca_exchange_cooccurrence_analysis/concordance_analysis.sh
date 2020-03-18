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
source ${PROGRAM_DIR}/vcfvenv3.8/bin/activate

#####################################
## FORMAT RAW SAMPLE GENOTYPE VCFS ##
echo "formatting raw sample genotype vcfs"
cd ${WORK_DIR}
# Extract brca1 region (GRCh38 coordinates)
if [ ! -f "${WORK_DIR}/raw_sample.phased.brca1region.vcf.gz" ]; then
    ${PROGRAM_DIR}/bcftools filter --threads 8 --regions chr17:43044295-43170245 ${CHR17_BCF} > raw_sample.phased.brca1region.vcf
    bgzip raw_sample.phased.brca1region.vcf
    tabix -p vcf raw_sample.phased.brca1region.vcf.gz
fi
# Extract brca2 region (GRCh38 coordinates)
if [ ! -f "${WORK_DIR}/raw_sample.phased.brca2region.vcf.gz" ]; then
    ${PROGRAM_DIR}/bcftools filter --threads 8 --regions chr13:32310000-32400000 ${CHR13_BCF} > raw_sample.phased.brca2region.vcf
    bgzip raw_sample.phased.brca2region.vcf
    tabix -p vcf raw_sample.phased.brca2region.vcf.gz
fi
echo "formatting raw sample genotype vcfs DONE"
echo ""

################################
## RUN VCF COMPARISON FILTERS ##
echo "running vcf normalization and comparison filters"
export HGVS_SEQREPO_DIR=${WORK_DIR}/seqrepo/2019-06-20
export UTA_DB_URL=postgresql://anonymous@localhost:15032/uta/uta_20170117
run_vcf_comparison () {
    BASE_VCF="$1"
    QUERY_VCF="$2"
    OUTPUT_FILENAME="$3"
    PROGRAM_DIR="$4"
    #python3.8 ${PROGRAM_DIR}/hgvs_normalize.py \
    #    -i ${BASE_VCF} \
    #    -o ${OUTPUT_FILENAME}.base_vcf.norm.vcf \
    #    -r hg38.p12.fa \
    #    -g ncbiRefSeq.txt
    #${PROGRAM_DIR}/vcf-sort -p 8 ${OUTPUT_FILENAME}.base_vcf.norm.vcf > ${OUTPUT_FILENAME}.base_vcf.norm.sorted.vcf
    #bgzip ${OUTPUT_FILENAME}.base_vcf.norm.sorted.vcf
    #tabix -p vcf ${OUTPUT_FILENAME}.base_vcf.norm.sorted.vcf.gz
    #rm -f ${OUTPUT_FILENAME}.base_vcf.norm.vcf
    
    ## Split VCFs into chunks and process each chunk in parallel
    #INPUT_VCF=${BASE_VCF}
    #BASENAME_VCF=$(basename "$INPUT_VCF")
    #VCF_FILENAME=${BASENAME_VCF%.vcf.gz}
    #NUM_CHUNKS=8
    #zcat ${INPUT_VCF} | tee >(head -n 10000 - | grep "^#" >${VCF_FILENAME}.header) | grep -v "^#" - >${VCF_FILENAME}.variants
    ##split into ${NUM_CHUNKS} chunks
    #split -d -n l/${NUM_CHUNKS} ${VCF_FILENAME}.variants ${VCF_FILENAME}
    ##reattach the header to each and clean up
    #for i in ${VCF_FILENAME}[0-9][0-9];do cat ${VCF_FILENAME}.header $i >$i.vcf && rm -f $i;done
    #rm -f ${VCF_FILENAME}.header ${VCF_FILENAME}.variants
    #normalize_base_vcf_pids=()
    #for i in ${NUM_CHUNKS}; do
    #    file_index=$((${i}-1))
    #    python3.8 ${PROGRAM_DIR}/hgvs_normalize.py \
    #        -i "${VCF_FILENAME}0${file_index}.vcf" \
    #        -o "${OUTPUT_FILENAME}.base_vcf.norm.${file_index}.vcf" \
    #        -r hg38.p12.fa \
    #        -g ncbiRefSeq.txt &
    #    normalize_base_vcf_pids[${i}]=$!
    #done
    #for pid in ${normalize_base_vcf_pids[*]}; do
    #    wait $pid
    #done
    #sorted_base_vcfs=()
    #for i in ${NUM_CHUNKS}; do
    #    file_index=$((${i}-1))
    #    ${PROGRAM_DIR}/vcf-sort \
    #        ${OUTPUT_FILENAME}.base_vcf.norm.${file_index}.vcf \
    #        > ${OUTPUT_FILENAME}.base_vcf.norm.sorted.${file_index}.vcf && rm -f ${OUTPUT_FILENAME}.base_vcf.norm.${file_index}.vcf
    #    sorted_base_vcfs+=(${OUTPUT_FILENAME}.base_vcf.norm.sorted.${file_index}.vcf)
    #done
    #${PROGRAM_DIR}/bcftools concat -O b ${sorted_base_vcfs[@]} -o ${OUTPUT_FILENAME}.base_vcf.norm.vcf.gz
    #tabix -p vcf ${OUTPUT_FILENAME}.base_vcf.norm.vcf.gz
    #rm ${sorted_base_vcfs[@]}
    
    #echo "Normalizing ${QUERY_VCF}"
    #if [ ! -f "${WORK_DIR}/${OUTPUT_FILENAME}.sorted.vcf.gz" ]; then
    #    # Separate vcf variant data from genotype data
    #    ${PROGRAM_DIR}/bcftools view -G ${QUERY_VCF} -O v -o ${OUTPUT_FILENAME}.query_vcf.no_genotypes.vcf
    #    ${PROGRAM_DIR}/bcftools query -f '[\t%GT]\n' ${QUERY_VCF} > ${OUTPUT_FILENAME}.query_vcf.genotypes.vcf
    #    bgzip ${OUTPUT_FILENAME}.query_vcf.no_genotypes.vcf
    #    tabix -p vcf ${OUTPUT_FILENAME}.query_vcf.no_genotypes.vcf.gz
    #    # Normalize only the variant data
    #    python3.8 ${PROGRAM_DIR}/hgvs_normalize.py \
    #        -i ${OUTPUT_FILENAME}.query_vcf.no_genotypes.vcf.gz \
    #        -o ${OUTPUT_FILENAME}.query_vcf.norm.no_genotypes.vcf \
    #        -r hg38.p12.fa \
    #        -g ncbiRefSeq.txt
    #    
    #    # Merge normalized variant data with genotype data
    #    ${PROGRAM_DIR}/bcftools view -H ${OUTPUT_FILENAME}.query_vcf.norm.no_genotypes.vcf > ${OUTPUT_FILENAME}.query_vcf.norm.no_genotypes.no_header.vcf
    #    paste ${OUTPUT_FILENAME}.query_vcf.norm.no_genotypes.no_header.vcf ${OUTPUT_FILENAME}.query_vcf.genotypes.vcf > ${OUTPUT_FILENAME}.query_vcf.norm.paste.vcf
    #    ${PROGRAM_DIR}/bcftools view -h ${OUTPUT_FILENAME}.query_vcf.norm.no_genotypes.vcf >> ${OUTPUT_FILENAME}.query_vcf.norm.vcf
    #    cat ${OUTPUT_FILENAME}.query_vcf.norm.paste.vcf >> ${OUTPUT_FILENAME}.query_vcf.norm.vcf
    #    ${PROGRAM_DIR}/vcf-sort -p 8 ${OUTPUT_FILENAME}.query_vcf.norm.vcf > ${OUTPUT_FILENAME}.query_vcf.norm.sorted.vcf
    #    bgzip ${OUTPUT_FILENAME}.query_vcf.norm.sorted.vcf
    #    tabix -p vcf ${OUTPUT_FILENAME}.query_vcf.norm.sorted.vcf.gz
    #    # Cleanup intermediate files
    #    rm -f  ${OUTPUT_FILENAME}.query_vcf.norm.vcf ${OUTPUT_FILENAME}.query_vcf.norm.paste.vcf ${OUTPUT_FILENAME}.query_vcf.norm.no_genotypes.vcf ${OUTPUT_FILENAME}.query_vcf.norm.no_genotypes.no_header.vcf ${OUTPUT_FILENAME}.query_vcf.genotypes.vcf ${OUTPUT_FILENAME}.query_vcf.no_genotypes.vcf.gz ${OUTPUT_FILENAME}.query_vcf.no_genotypes.vcf.gz.tbi
    #fi
    #INPUT_VCF=${QUERY_VCF}
    #BASENAME_VCF=$(basename "$INPUT_VCF")
    #VCF_FILENAME=${BASENAME_VCF%.vcf.gz}
    #NUM_CHUNKS=8
    #zcat ${INPUT_VCF} | tee >(head -n 10000 - | grep "^#" >${VCF_FILENAME}.header) | grep -v "^#" - >${VCF_FILENAME}.variants
    ##split into ${NUM_CHUNKS} chunks
    #split -d -n l/${NUM_CHUNKS} ${VCF_FILENAME}.variants ${VCF_FILENAME}
    ##reattach the header to each and clean up
    #for i in ${VCF_FILENAME}[0-9][0-9];do cat ${VCF_FILENAME}.header $i >$i.vcf && rm -f $i;done
    #rm -f ${VCF_FILENAME}.header ${VCF_FILENAME}.variants
    #normalize_query_vcf_pids=()
    #for i in ${NUM_CHUNKS}; do
    #    file_index=$((${i}-1))
    #    python3.8 ${PROGRAM_DIR}/hgvs_normalize.py \
    #        -i "${VCF_FILENAME}0${file_index}.vcf" \
    #        -o "${OUTPUT_FILENAME}.query_vcf.norm.${file_index}.vcf" \
    #        -r hg38.p12.fa \
    #        -g ncbiRefSeq.txt &
    #    normalize_query_vcf_pids[${i}]=$!
    #done
    #for pid in ${normalize_query_vcf_pids[*]}; do
    #    wait $pid
    #done
    #sorted_query_vcfs=()
    #for i in ${NUM_CHUNKS}; do
    #    file_index=$((${i}-1))
    #    ${PROGRAM_DIR}/vcf-sort \
    #        ${OUTPUT_FILENAME}.query_vcf.norm.${file_index}.vcf \
    #        > ${OUTPUT_FILENAME}.query_vcf.norm.sorted.${file_index}.vcf && rm -f ${OUTPUT_FILENAME}.query_vcf.norm.${file_index}.vcf
    #    sorted_query_vcfs+=(${OUTPUT_FILENAME}.query_vcf.norm.sorted.${file_index}.vcf)
    #done
    #${PROGRAM_DIR}/bcftools concat -O b ${sorted_query_vcfs[@]} -o ${OUTPUT_FILENAME}.query_vcf.norm.vcf.gz
    #tabix -p vcf ${OUTPUT_FILENAME}.query_vcf.norm.vcf.gz
    #rm ${sorted_query_vcfs[@]}
    
    ${PROGRAM_DIR}/bcftools isec \
        -O v \
        -n =2 -w 1 \
        ${QUERY_VCF} \
        ${BASE_VCF} \
        > ${OUTPUT_FILENAME}.vcf
    ${PROGRAM_DIR}/vcf-sort -p 8 ${OUTPUT_FILENAME}.vcf > ${OUTPUT_FILENAME}.sorted.vcf
    bgzip ${OUTPUT_FILENAME}.sorted.vcf
    tabix -p vcf ${OUTPUT_FILENAME}.sorted.vcf.gz
    rm -f ${OUTPUT_FILENAME}.vcf
}

run_normalize_vcf () {
    INPUT_VCF="$1"
    OUTPUT_FILENAME="$2"
    SAMPLE_GENOTYPES="$3"
    PROGRAM_DIR="$4"
    echo "Normalizing ${INPUT_VCF}"
    if [ ${SAMPLE_GENOTYPES} = 'false' ]; then
        python3.8 ${PROGRAM_DIR}/hgvs_normalize.py \
            -i ${INPUT_VCF} \
            -o ${OUTPUT_FILENAME}.norm.vcf \
            -r hg38.p12.fa \
            -g ncbiRefSeq.txt
        ${PROGRAM_DIR}/vcf-sort -p 8 ${OUTPUT_FILENAME}.norm.vcf > ${OUTPUT_FILENAME}.norm.sorted.vcf
        bgzip ${OUTPUT_FILENAME}.norm.sorted.vcf
        tabix -p vcf ${OUTPUT_FILENAME}.norm.sorted.vcf.gz
        rm -f ${OUTPUT_FILENAME}.norm.vcf
    else
        # Separate vcf variant data from genotype data
        ${PROGRAM_DIR}/bcftools view -G ${INPUT_VCF} -O v -o ${OUTPUT_FILENAME}.no_genotypes.vcf
        ${PROGRAM_DIR}/bcftools query -f '[\t%GT]\n' ${INPUT_VCF} > ${OUTPUT_FILENAME}.genotypes.vcf
        bgzip ${OUTPUT_FILENAME}.no_genotypes.vcf
        tabix -p vcf ${OUTPUT_FILENAME}.no_genotypes.vcf.gz
        # Normalize only the variant data
        python3.8 ${PROGRAM_DIR}/hgvs_normalize.py \
            -i ${OUTPUT_FILENAME}.no_genotypes.vcf.gz \
            -o ${OUTPUT_FILENAME}.norm.no_genotypes.vcf \
            -r hg38.p12.fa \
            -g ncbiRefSeq.txt
        
        # Merge normalized variant data with genotype data
        ${PROGRAM_DIR}/bcftools view -H ${OUTPUT_FILENAME}.norm.no_genotypes.vcf > ${OUTPUT_FILENAME}.norm.no_genotypes.no_header.vcf
        paste ${OUTPUT_FILENAME}.norm.no_genotypes.no_header.vcf ${OUTPUT_FILENAME}.genotypes.vcf > ${OUTPUT_FILENAME}.norm.paste.vcf
        ${PROGRAM_DIR}/bcftools view -h ${OUTPUT_FILENAME}.norm.no_genotypes.vcf >> ${OUTPUT_FILENAME}.norm.vcf
        cat ${OUTPUT_FILENAME}.norm.paste.vcf >> ${OUTPUT_FILENAME}.norm.vcf
        ${PROGRAM_DIR}/vcf-sort -p 8 ${OUTPUT_FILENAME}.norm.vcf > ${OUTPUT_FILENAME}.norm.sorted.vcf
        bgzip ${OUTPUT_FILENAME}.norm.sorted.vcf
        tabix -p vcf ${OUTPUT_FILENAME}.norm.sorted.vcf.gz
        # Cleanup intermediate files
        rm -f  ${OUTPUT_FILENAME}.norm.vcf ${OUTPUT_FILENAME}.norm.paste.vcf ${OUTPUT_FILENAME}.norm.no_genotypes.vcf ${OUTPUT_FILENAME}.norm.no_genotypes.no_header.vcf ${OUTPUT_FILENAME}.genotypes.vcf ${OUTPUT_FILENAME}.no_genotypes.vcf.gz ${OUTPUT_FILENAME}.no_genotypes.vcf.gz.tbi
    fi
}

# Normalize VCFs
INPUT_VCF="${WORK_DIR}/brcaexchange_variants.brca1.pathogenic.vcf.gz"
OUTPUT_FILENAME="brca1_exchange_PATHOGENIC"
SAMPLE_GENOTYPES="false"
run_normalize_vcf ${INPUT_VCF} ${OUTPUT_FILENAME} ${SAMPLE_GENOTYPES} ${PROGRAM_DIR}

INPUT_VCF="${WORK_DIR}/brcaexchange_variants.brca2.pathogenic.vcf.gz"
OUTPUT_FILENAME="brca2_exchange_PATHOGENIC"
SAMPLE_GENOTYPES="false"
run_normalize_vcf ${INPUT_VCF} ${OUTPUT_FILENAME} ${SAMPLE_GENOTYPES} ${PROGRAM_DIR}

INPUT_VCF="${WORK_DIR}/brcaexchange_variants.brca1.vus.vcf.gz"
OUTPUT_FILENAME="brca1_exchange_VUS"
SAMPLE_GENOTYPES="false"
run_normalize_vcf ${INPUT_VCF} ${OUTPUT_FILENAME} ${SAMPLE_GENOTYPES} ${PROGRAM_DIR}

INPUT_VCF="${WORK_DIR}/brcaexchange_variants.brca2.vus.vcf.gz"
OUTPUT_FILENAME="brca2_exchange_VUS"
SAMPLE_GENOTYPES="false"
run_normalize_vcf ${INPUT_VCF} ${OUTPUT_FILENAME} ${SAMPLE_GENOTYPES} ${PROGRAM_DIR}

INPUT_VCF="${WORK_DIR}/raw_sample.phased.brca1region.vcf.gz"
OUTPUT_FILENAME="brca1_sample_genotypes"
SAMPLE_GENOTYPES="true"
run_normalize_vcf ${INPUT_VCF} ${OUTPUT_FILENAME} ${SAMPLE_GENOTYPES} ${PROGRAM_DIR}

INPUT_VCF="${WORK_DIR}/raw_sample.phased.brca2region.vcf.gz"
OUTPUT_FILENAME="brca2_sample_genotypes"
SAMPLE_GENOTYPES="true"
run_normalize_vcf ${INPUT_VCF} ${OUTPUT_FILENAME} ${SAMPLE_GENOTYPES} ${PROGRAM_DIR}

# Process brca1 pathogenic variants
BASE_VCF="${WORK_DIR}/brca1_exchange_PATHOGENIC.norm.sorted.vcf.gz"
QUERY_VCF="${WORK_DIR}/brca1_sample_genotypes.norm.sorted.vcf.gz"
OUTPUT_FILENAME="brca1_PATHOGENIC"
run_vcf_comparison ${BASE_VCF} ${QUERY_VCF} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

# Process brca1 vus variants
BASE_VCF="${WORK_DIR}/brca1_exchange_VUS.norm.sorted.vcf.gz"
QUERY_VCF="${WORK_DIR}/brca1_sample_genotypes.norm.sorted.vcf.gz"
OUTPUT_FILENAME="brca1_VUS"
run_vcf_comparison ${BASE_VCF} ${QUERY_VCF} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

# Process brca2 pathogenic variants
BASE_VCF="${WORK_DIR}/brca2_exchange_PATHOGENIC.norm.sorted.vcf.gz"
QUERY_VCF="${WORK_DIR}/brca2_sample_genotypes.norm.sorted.vcf.gz"
OUTPUT_FILENAME="brca2_PATHOGENIC"
run_vcf_comparison ${BASE_VCF} ${QUERY_VCF} ${OUTPUT_FILENAME} ${PROGRAM_DIR}

# Process brca2 vus variants
BASE_VCF="${WORK_DIR}/brca2_exchange_VUS.norm.sorted.vcf.gz"
QUERY_VCF="${WORK_DIR}/brca2_sample_genotypes.norm.sorted.vcf.gz"
OUTPUT_FILENAME="brca2_VUS"
run_vcf_comparison ${BASE_VCF} ${QUERY_VCF} ${OUTPUT_FILENAME} ${PROGRAM_DIR}
echo "running vcf comparison filters DONE"
echo ""


##################################
## RUN THE CONCORDANCE ANALYSIS ##
echo "running the concordance analysis"
python3.8 ${PROGRAM_DIR}/detect_vus_benign.py -i "${WORK_DIR}/brca1_VUS.sorted.vcf.gz" -j "${WORK_DIR}/brca1_PATHOGENIC.sorted.vcf.gz" -o "${WORK_DIR}/brca1_coocurrence_report.txt" -v "${WORK_DIR}/brca1_apparent_benign_VUS.vcf"
python3.8 ${PROGRAM_DIR}/detect_vus_benign.py -i "${WORK_DIR}/brca2_VUS.sorted.vcf.gz" -j "${WORK_DIR}/brca2_PATHOGENIC.sorted.vcf.gz" -o "${WORK_DIR}/brca2_coocurrence_report.txt" -v "${WORK_DIR}/brca2_apparent_benign_VUS.vcf"

${PROGRAM_DIR}/vcf-sort ${WORK_DIR}/brca1_apparent_benign_VUS.vcf > ${WORK_DIR}/brca1_apparent_benign_VUS.sorted.vcf
bgzip ${WORK_DIR}/brca1_apparent_benign_VUS.sorted.vcf
tabix -p vcf ${WORK_DIR}/brca1_apparent_benign_VUS.sorted.vcf.gz
rm ${WORK_DIR}/brca1_apparent_benign_VUS.vcf
${PROGRAM_DIR}/vcf-sort ${WORK_DIR}/brca2_apparent_benign_VUS.vcf > ${WORK_DIR}/brca2_apparent_benign_VUS.sorted.vcf
bgzip ${WORK_DIR}/brca2_apparent_benign_VUS.sorted.vcf
tabix -p vcf ${WORK_DIR}/brca2_apparent_benign_VUS.sorted.vcf.gz
rm ${WORK_DIR}/brca2_apparent_benign_VUS.vcf
deactivate
echo "running the concordance analysis DONE"
echo ""

exit

