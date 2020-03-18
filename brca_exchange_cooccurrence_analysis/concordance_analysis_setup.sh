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
if [ $# -lt 2 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## Parse through arguments
while getopts "p:w:h" OPTION; do
    case $OPTION in
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

# Install a new version of sqlite3 (for seqrepo compatibility)
cd ${WORK_DIR}
wget https://www.sqlite.org/2020/sqlite-autoconf-3310100.tar.gz
tar xvfz sqlite-autoconf-3310100.tar.gz
cd sqlite-autoconf-3310100
mkdir -p ${WORK_DIR}/lib
mkdir -p ${WORK_DIR}/bin
./configure --prefix=${WORK_DIR}
make
make install

# Install a new version of python3
cd ${WORK_DIR}
wget https://www.python.org/ftp/python/3.8.2/Python-3.8.2.tgz
tar xvfz Python-3.8.2.tgz
cd Python-3.8.2
LD_RUN_PATH=${WORK_DIR}/lib ./configure --enable-loadable-sqlite-extensions --prefix=${WORK_DIR} --with-pydebug
perl -pi.orig -e "s|(?<=sqlite_inc_paths = )\[|['${WORK_DIR}/include',\n|" setup.py
LD_RUN_PATH=${WORK_DIR}/lib make
LD_RUN_PATH=${WORK_DIR}/lib make altinstall

cd ${WORK_DIR}
${WORK_DIR}/bin/python3.8 -m venv ${PROGRAM_DIR}/vcfvenv3.8
source ${PROGRAM_DIR}/vcfvenv3.8/bin/activate
python3.8 -m pip install --upgrade pip
python3.8 -m pip install cython
python3.8 -m pip install IPython==5.0
python3.8 -m pip install hgvs
python3.8 -m pip install requests pyvcf hgvs pyhgvs pyfaidx
python3.8 -m pip install biocommons.seqrepo
echo "building python virtualenv DONE"
echo ""

##############################################
## download brcaexchange data to vcf format ##
echo "downloading and formatting brcaexchange data"
cd ${WORK_DIR}
python3.8 ${PROGRAM_DIR}/extract_brcaexchange_data.py -o brcaexchange_variants
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
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.fa.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
gzip -df hg38.p12.fa.gz > hg38.p12.fa
gzip -df ncbiRefSeq.txt.gz > ncbiRefSeq.txt
faidx --no-output hg38.p12.fa
if [ $(docker container ps --filter name="uta_20170117" | wc -l) == 1 ]; then
    docker run -d --name uta_20170117 -p 15032:5432 biocommons/uta:uta_20170117
fi
mkdir ${WORK_DIR}/seqrepo
seqrepo -r ${WORK_DIR}/seqrepo pull -i 2019-06-20
echo "formatting reference index DONE"
echo ""
deactivate

exit

