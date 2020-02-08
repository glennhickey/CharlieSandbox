## Running instructions

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

