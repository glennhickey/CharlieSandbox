## Running instructions

### Running the mosaicism snp extraction script

cd into the home directory and download the code.
```
cd ~
git clone https://github.com/cmarkello/Sandbox.git
```

Create a python virtual environment and install pyvcf
```
python3 -m venv pyvcf_venv
source pyvcf_venv/bin/activate
pip install pyvcf
deactivate
```

Launch an interactive node and setup input variables
```
sinteractive
module load python/3.7
source ~/pyvcf_venv/bin/activate
INPUT_JOINT_VCF='/path/to/cohort.vcf.gz'
MATERNAL_SAMPLE_NAME='HG004'
PATERNAL_SAMPLE_NAME='HG003'
OUTPUT_BASENAME='mosaicism_list'
```

Run the script
```
time python3 ~/Sandbox/mosaicism_analysis/extract_good_mosaic_variants.py -i ${INPUT_JOINT_VCF} -m ${MATERNAL_SAMPLE_NAME} -p ${PATERNAL_SAMPLE_NAME} -o ${OUTPUT_BASENAME}
```

Output will be found in the file `${SIBLING_NAME}.${OUTPUT_BASENAME}.tsv` for each sibling id in the cohort.
Each file is in tab delimited format where the columns represent in-order: Chromosome, Position, maternal allele depth for sibling, paternal allele depth for sibling

### Running the mosaicism snp extraction script for a specific cohort in a population vcf

To run the script on a joint-called vcf that contains more than one cohort, you will want to specify the sibling ids exlicitly as they are represented in the vcf by using the `-s` parameter for each sibling.
```
INPUT_JOINT_VCF='/path/to/cohort.vcf.gz'
MATERNAL_SAMPLE_NAME='HG004'
PATERNAL_SAMPLE_NAME='HG003'
OUTPUT_BASENAME='mosaicism_list'
SIBLING_1_NAME='HG002'
SIBLING_2_NAME='NA12878'

time python3 ~/Sandbox/mosaicism_analysis/extract_good_mosaic_variants.py -i ${INPUT_JOINT_VCF} -m ${MATERNAL_SAMPLE_NAME} -p ${PATERNAL_SAMPLE_NAME} -s ${SIBLING_1_NAME} -s ${SIBLING_2_NAME} -o ${OUTPUT_BASENAME}
```
