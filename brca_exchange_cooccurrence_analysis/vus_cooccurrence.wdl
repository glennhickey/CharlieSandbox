version 1.0

workflow vus_cooccurrence {
    input {
        File SAMPLE_BCF
        String OUTPUT_NAME
        String REGION
    }
    call setup_brcaexchange_data {
        in_region=REGION,
        outname=OUTPUT_NAME
    }
    output {
    }
}

task setup_brcaexchange_data {
    input {
        String in_region
        String outname
    }
    command <<<
        set -exu -o pipefail
        
        wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.fa.gz
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
        
        IFS=':' read -ra split_region <<< "~{in_region}"
        chromosome_name="${split_region[0]}"
        out_brcaexchange_path_vcf=""
        out_brcaexchange_vus_vcf=""
        if [ ${chromosome_name} = "chr17" ]; then
            python3 ./extract_brcaexchange_data.py -g BRCA1 -o ~{outname}
            out_brcaexchange_path_vcf="~{outname}.brca1.pathogenic.vcf"
            out_brcaexchange_uvs_vcf="~{outname}.brca1.vus.vcf"
        elif [ ${chromosome_name} = "chr13" ]; then
            python3 ./extract_brcaexchange_data.py -g BRCA2 -o ~{outname}
            out_brcaexchange_path_vcf="~{outname}.brca2.pathogenic.vcf"
            out_brcaexchange_vus_vcf="~{outname}.brca2.vus.vcf"
        fi
    >>>
    output {
        File reference_file = "hg38.p12.fa.gz"
        File ref_seq_file = "ncbiRefSeq.txt.gz"
        File pathogenic_vcf = glob("*.pathogenic.vcf")[0]
        File vus_vcf = glob("*.vus.vcf")[0]
    }
    runtime {
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

