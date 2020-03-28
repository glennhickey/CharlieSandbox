## Running instructions

### Running the workflow using the Docker container and Cromwell + WDL

cd into the home directory and download the code.
```
cd ~
git clone https://github.com/cmarkello/Sandbox.git
cd Sandbox/brca_exchange_cooccurrence_analysis
```
Create a python virtual environment and install miniwdl
```
python3 -m venv miniwdl_venv
source miniwdl_venv/bin/activate
pip install miniwdl
```
Setup input variables to the full paths to the population joint-genotyped and phased vcfs along with the analysis region in chromosomal coordinates. The following are example input.
```
INPUT_SAMPLE_BCF='/path/to/chr17.phased.bcf'
INPUT_SAMPLE_BCF_INDEX='/path/to/chr17.phased.bcf.csi'
OUTPUT_NAME='test_brca1_wdl'
REGION='chr17:43044295-43170245'
```
Run the WDL workflow on Cromwell using miniwdl
```
miniwdl cromwell ./Sandbox/brca_exchange_cooccurrence_analysis/vus_cooccurrence.wdl \
SAMPLE_BCF=${INPUT_SAMPLE_BCF} \
SAMPLE_BCF_INDEX=${INPUT_SAMPLE_BCF_INDEX} \
OUTPUT_NAME=${OUTPUT_NAME} \
REGION=${REGION} \
-c ./Sandbox/brca_exchange_cooccurrence_analysis/cromwell.local.conf \
-d '/path/to/test_final_outputs'
```
The workflow with output the following:
- Apparent Benign VUS variants VCF located in `/path/to/test_final_outputs/*_vus_cooccurrence/output_links/apparent_benign_vus_vcf`
- Basic cooccurrence count report located in `/path/to/test_final_outputs/*_vus_cooccurrence/output_links/cooccurrence_report`

### Running the workflow using scripts directly

cd into the home directory and download the code.
```
cd ~
git clone https://github.com/cmarkello/Sandbox.git
cd Sandbox/brca_exchange_cooccurrence_analysis
```

Requisite programs and their respective sources
```
extract_brcaexchange_data.py    (from https://github.com/cmarkello/Sandbox/tree/master/brca_exchange_cooccurrence_analysis)
detect_vus_benign.py            (from https://github.com/cmarkello/Sandbox/tree/master/brca_exchange_cooccurrence_analysis)
bcftools                        (from http://www.htslib.org/download/)
rtg                             (from https://github.com/RealTimeGenomics/rtg-tools)
vcf-sort                        (from https://vcftools.github.io/perl_module.html)
```

To run the program simply give it the full paths to the population joint-genotyped and phased vcfs along with its requisite programs
and working directory.
```
CHR17_BCF="/path/to/chr17.phased.bcf"
CHR13_BCF="/path/to/chr13.phased.bcf"
PROGRAM_DIR="/path/to/programs/directory"
WORK_DIR="/path/to/workflow/where/program/and/data/are/exicuted"
```

Run the program
```
time ${PROGRAM_DIR}/concordance_analysis_setup.sh -a ${CHR17_BCF} -b ${CHR13_BCF} -p ${PROGRAM_DIR} -w ${WORK_DIR}
```

