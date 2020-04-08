version 1.0

workflow vus_cooccurrence {
    input {
        File SAMPLE_BCF
        File SAMPLE_BCF_INDEX
        String OUTPUT_NAME
        String REGION
    }
    
    call setup_brcaexchange_data {
        input:
            in_region=REGION,
            outname=OUTPUT_NAME
    }
    call normalize_brcaexchange_data as normalize_brca_path {
        input:
            in_vcf=setup_brcaexchange_data.pathogenic_vcf,
            in_ref_file=setup_brcaexchange_data.reference_file,
            in_seq_file=setup_brcaexchange_data.ref_seq_file
    }
    call normalize_brcaexchange_data as normalize_brca_vus {
        input:
            in_vcf=setup_brcaexchange_data.vus_vcf,
            in_ref_file=setup_brcaexchange_data.reference_file,
            in_seq_file=setup_brcaexchange_data.ref_seq_file
    }
    
    call setup_sample_data {
        input:
            in_sample_bcf=SAMPLE_BCF,
            in_sample_bcf_index=SAMPLE_BCF_INDEX,
            in_region=REGION,
            in_ref_file=setup_brcaexchange_data.reference_file,
            in_seq_file=setup_brcaexchange_data.ref_seq_file
    }
    
    call intersect_variants as intersect_path_variants {
        input:
            in_base_vcf=normalize_brca_path.normalized_vcf,
            in_base_vcf_index=normalize_brca_path.normalized_vcf_index,
            in_query_vcf=setup_sample_data.normalized_vcf,
            in_query_vcf_index=setup_sample_data.normalized_vcf_index
    }
    call intersect_variants as intersect_vus_variants {
        input:
            in_base_vcf=normalize_brca_vus.normalized_vcf,
            in_base_vcf_index=normalize_brca_vus.normalized_vcf_index,
            in_query_vcf=setup_sample_data.normalized_vcf,
            in_query_vcf_index=setup_sample_data.normalized_vcf_index
    }
    
    call dectect_vus_benign {
        input:
            in_intersect_vus_vcf=intersect_vus_variants.intersected_vcf,
            in_intersect_vus_vcf_index=intersect_vus_variants.intersected_vcf_index,
            in_intersect_path_vcf=intersect_path_variants.intersected_vcf,
            in_intersect_path_vcf_index=intersect_path_variants.intersected_vcf_index,
            outname=OUTPUT_NAME
    }
    
    output {
        File cooccurrence_report = dectect_vus_benign.cooccurrence_report
        File complete_cooccurrence_report = dectect_vus_benign.complete_cooccurrence_report
        File apparent_benign_vus_vcf = dectect_vus_benign.apparent_benign_vus_vcf
        File apparent_benign_vus_vcf_index = dectect_vus_benign.apparent_benign_vus_vcf_index
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
        gzip -d hg38.p12.fa.gz
        gzip -d ncbiRefSeq.txt.gz
        
        IFS=':' read -ra split_region <<< "~{in_region}"
        chromosome_name="${split_region[0]}"
        out_brcaexchange_path_vcf=""
        out_brcaexchange_vus_vcf=""
        ls -l /usr/src/app/
        if [ ${chromosome_name} = "chr17" ]; then
            python3 /usr/src/app/extract_brcaexchange_data.py -g BRCA1 -o ~{outname}
            out_brcaexchange_path_vcf="~{outname}.brca1.pathogenic.vcf"
            out_brcaexchange_uvs_vcf="~{outname}.brca1.vus.vcf"
        elif [ ${chromosome_name} = "chr13" ]; then
            python3 /usr/src/app/extract_brcaexchange_data.py -g BRCA2 -o ~{outname}
            out_brcaexchange_path_vcf="~{outname}.brca2.pathogenic.vcf"
            out_brcaexchange_vus_vcf="~{outname}.brca2.vus.vcf"
        fi
    >>>
    output {
        File reference_file = "hg38.p12.fa"
        File ref_seq_file = "ncbiRefSeq.txt"
        File pathogenic_vcf = glob("*.pathogenic.vcf")[0]
        File vus_vcf = glob("*.vus.vcf")[0]
    }
    runtime {
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

task normalize_brcaexchange_data {
    input {
        File in_vcf
        File in_ref_file
        File in_seq_file
    }
    command <<<
        set -exu -o pipefail
        
        ln -s ~{in_vcf} input.vcf
        bgzip input.vcf
        tabix -p vcf input.vcf.gz
       
        mkdir -p ${PWD}/seqrepo
        seqrepo -r ${PWD}/seqrepo pull -i 2019-06-20
        export HGVS_SEQREPO_DIR=${PWD}/seqrepo/2019-06-20
        python3 /usr/src/app/hgvs_normalize.py \
            -i input.vcf.gz \
            -o brcaexchange.norm.vcf \
            -r ~{in_ref_file} \
            -g ~{in_seq_file}
        
        vcf-sort -p 8 brcaexchange.norm.vcf > brcaexchange.norm.sorted.vcf
        bgzip brcaexchange.norm.sorted.vcf
        tabix -p vcf brcaexchange.norm.sorted.vcf.gz
        rm -f brcaexchange.norm.vcf
    >>>
    output {
        File normalized_vcf = "brcaexchange.norm.sorted.vcf.gz"
        File normalized_vcf_index = "brcaexchange.norm.sorted.vcf.gz.tbi"
    }
    runtime {
        cpu: 8
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

task setup_sample_data {
    input {
        File in_sample_bcf
        File in_sample_bcf_index
        String in_region
        File in_ref_file
        File in_seq_file
    }
    command <<<
        set -exu -o pipefail
        
        # Extract specific region from vcf
        ln -s ~{in_sample_bcf} input_raw_sample.bcf
        ln -s ~{in_sample_bcf_index} input_raw_sample.bcf.csi
        bcftools filter --threads 8 --regions ~{in_region} input_raw_sample.bcf > raw_sample.vcf
        bgzip raw_sample.vcf
        tabix -p vcf raw_sample.vcf.gz
        
        # Separate vcf variant data from genotype data
        bcftools view -G raw_sample.vcf.gz -O v -o raw_sample.no_genotypes.vcf
        bcftools query -f 'GT\t[%GT\t]\n' raw_sample.vcf.gz > raw_sample.genotypes.vcf
        bgzip raw_sample.no_genotypes.vcf
        tabix -p vcf raw_sample.no_genotypes.vcf.gz
        
        # Normalize only the variant data
        mkdir -p ${PWD}/seqrepo
        seqrepo -r ${PWD}/seqrepo pull -i 2019-06-20
        export HGVS_SEQREPO_DIR=${PWD}/seqrepo/2019-06-20
        python3 /usr/src/app/hgvs_normalize.py \
            -i raw_sample.no_genotypes.vcf.gz \
            -o raw_sample.no_genotypes.norm.vcf \
            -r ~{in_ref_file} \
            -g ~{in_seq_file}
        
        # Merge normalized variant data with genotype data
        bcftools view -H raw_sample.no_genotypes.norm.vcf > raw_sample.no_genotypes.norm.no_header.vcf
        paste raw_sample.no_genotypes.norm.no_header.vcf raw_sample.genotypes.vcf > raw_sample.norm.paste.vcf
        bcftools view -h raw_sample.vcf.gz >> raw_sample.norm.vcf
        cat raw_sample.norm.paste.vcf >> raw_sample.norm.vcf
        vcf-sort -p 8 raw_sample.norm.vcf > raw_sample.norm.sorted.vcf
        bgzip raw_sample.norm.sorted.vcf
        tabix -p vcf raw_sample.norm.sorted.vcf.gz
        rm -f raw_sample.norm.vcf raw_sample.norm.paste.vcf raw_sample.no_genotypes.norm.no_header.vcf raw_sample.no_genotypes.norm.vcf raw_sample.no_genotypes.vcf.gz raw_sample.no_genotypes.vcf.gz.tbi raw_sample.genotypes.vcf
    >>>
    output {
        File normalized_vcf = "raw_sample.norm.sorted.vcf.gz"
        File normalized_vcf_index = "raw_sample.norm.sorted.vcf.gz.tbi"
    }
    runtime {
        cpu: 8
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

task intersect_variants {
    input {
        File in_base_vcf
        File in_base_vcf_index
        File in_query_vcf
        File in_query_vcf_index
    }
    command <<<
        set -exu -o pipefail
        
        ln -s ~{in_base_vcf} base.vcf.gz
        ln -s ~{in_base_vcf_index} base.vcf.gz.tbi
        ln -s ~{in_query_vcf} query.vcf.gz
        ln -s ~{in_query_vcf_index} query.vcf.gz.tbi
        
        bcftools isec \
            -O v \
            -n =2 -w 1 \
            -o intersected_variants.vcf \
            query.vcf.gz \
            base.vcf.gz
        vcf-sort -p 8 intersected_variants.vcf > intersected_variants.sorted.vcf
        bgzip intersected_variants.sorted.vcf
        tabix -p vcf intersected_variants.sorted.vcf.gz
        #rm -f intersected_variants.vcf
    >>>
    output {
        File intersected_vcf = "intersected_variants.sorted.vcf.gz"
        File intersected_vcf_index = "intersected_variants.sorted.vcf.gz.tbi"
    }
    runtime {
        cpu: 8
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}

task dectect_vus_benign {
    input {
        File in_intersect_vus_vcf
        File in_intersect_vus_vcf_index
        File in_intersect_path_vcf
        File in_intersect_path_vcf_index
        String outname
    }
    command <<<
        set -exu -o pipefail
        
        ln -s ~{in_intersect_vus_vcf} vus.vcf.gz
        ln -s ~{in_intersect_vus_vcf_index} vus.vcf.gz.tbi
        ln -s ~{in_intersect_path_vcf} path.vcf.gz
        ln -s ~{in_intersect_path_vcf_index} path.vcf.gz.tbi
        
        python3 /usr/src/app/detect_vus_benign.py \
            -i vus.vcf.gz \
            -j path.vcf.gz \
            -o ~{outname}.cooccurrence_report.txt \
            -v ~{outname}.apperent_benign_vus_list.vcf
        
        vcf-sort -p 8 ~{outname}.apperent_benign_vus_list.vcf > ~{outname}.apperent_benign_vus_list.sorted.vcf
        bgzip ~{outname}.apperent_benign_vus_list.sorted.vcf
        tabix -p vcf ~{outname}.apperent_benign_vus_list.sorted.vcf.gz
        rm -f ~{outname}.apperent_benign_vus_list.vcf
    >>>
    output {
        File cooccurrence_report = "~{outname}.cooccurrence_report.txt"
        File complete_cooccurrence_report = "complete_~{outname}.cooccurrence_report.txt"
        File apparent_benign_vus_vcf = "~{outname}.apperent_benign_vus_list.sorted.vcf.gz"
        File apparent_benign_vus_vcf_index = "~{outname}.apperent_benign_vus_list.sorted.vcf.gz.tbi"
    }
    runtime {
        cpu: 8
        docker: 'quay.io/cmarkello/vus_cooccurrence:latest'
    }
}





