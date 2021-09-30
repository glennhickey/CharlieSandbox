
version 1.0

### deepvariant.wdl ###
## Author: Charles Markello
## Description: Core DeepVariant calling workflow for single sample datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow deepvariant {
    input {
        String SAMPLE_NAME                              # The sample name
        Array[String] CONTIGS                           # (OPTIONAL) Desired reference genome contigs, which are all paths in the XG index.
        File REFERENCE_FASTA_FILE                       # (OPTIONAL) Use this refernce instead of extracting paths from the XG. Required if the graph does not contain the entire reference (ex when making GRCh38 calls on a CHM13-based graph)
        File PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index. If not  given, paths are extracted from the XG and subset to chromosome-looking paths.
        File BAM_FILE
        File BAM_INDEX_FILE
        File? TRUTH_VCF                                 # Path to .vcf.gz to compare against
        File? TRUTH_VCF_INDEX                           # Path to Tabix index for TRUTH_VCF
        File? EVALUATION_REGIONS_BED                    # BED to restrict comparison against TRUTH_VCF to
        File? DV_MODEL_META                             # .meta file for a custom DeepVariant calling model
        File? DV_MODEL_INDEX                            # .index file for a custom DeepVariant calling model
        File? DV_MODEL_DATA                             # .data-00000-of-00001 file for a custom DeepVariant calling model
        Int MIN_MAPQ = 1                                # Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong.
        Int REALIGNMENT_EXPANSION_BASES = 160           # Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions.
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_DISK = 10
        Int MAP_CORES = 16
        Int MAP_DISK = 300
        Int MAP_MEM = 50
        Int CALL_CORES = 8
        Int CALL_DISK = 300
        Int CALL_MEM = 50
    }

    File reference_file = REFERENCE_FASTA_FILE
    
    call indexReference {
        input:
            in_reference_file=reference_file,
            in_index_disk=MAP_DISK,
            in_index_mem=MAP_MEM
    }
    File reference_index_file = indexReference.reference_index_file
    File reference_dict_file = indexReference.reference_dict_file

    # Split merged alignment by contigs list
    call splitBAMbyPath { 
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_bam_file=BAM_FILE,
            in_merged_bam_file_index=BAM_INDEX_FILE,
            in_path_list_file=PATH_LIST_FILE,
            in_map_cores=MAP_CORES,
            in_map_disk=MAP_DISK,
            in_map_mem=MAP_MEM
    }
 
    scatter (deepvariant_caller_input_files in zip(splitBAMbyPath.bam_contig_files, splitBAMbyPath.bam_contig_files_index)) {
        
        call runGATKRealignerTargetCreator {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=deepvariant_caller_input_files.left,
                in_bam_index_file=deepvariant_caller_input_files.right,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_reference_dict_file=reference_dict_file,
                in_call_disk=CALL_DISK
        }
        if (REALIGNMENT_EXPANSION_BASES != 0) {
            # We want the realignment targets to be wider
            call widenRealignmentTargets {
                input:
                    in_target_bed_file=runGATKRealignerTargetCreator.realigner_target_bed,
                    in_reference_index_file=reference_index_file,
                    in_expansion_bases=REALIGNMENT_EXPANSION_BASES,
                    in_call_disk=CALL_DISK
            }
        }
        File target_bed_file = select_first([widenRealignmentTargets.output_target_bed_file, runGATKRealignerTargetCreator.realigner_target_bed])
        call runAbraRealigner {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=deepvariant_caller_input_files.left,
                in_bam_index_file=deepvariant_caller_input_files.right,
                in_target_bed_file=target_bed_file,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_call_disk=CALL_DISK
        }
        call runDeepVariant {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=runAbraRealigner.indel_realigned_bam,
                in_bam_file_index=runAbraRealigner.indel_realigned_bam_index,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_model_meta_file=DV_MODEL_META,
                in_model_index_file=DV_MODEL_INDEX,
                in_model_data_file=DV_MODEL_DATA,
                in_min_mapq=MIN_MAPQ,
                in_call_cores=CALL_CORES,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
    }
    # Merge distributed variant called VCFs
    call concatClippedVCFChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_clipped_vcf_chunk_files=runDeepVariant.output_vcf_file,
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    # Extract either the normal or structural variant based VCFs and compress them
    call bgzipMergedVCF {
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_vcf_file=concatClippedVCFChunks.output_merged_vcf,
            in_vg_container="quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4",
            in_call_disk=CALL_DISK,
            in_call_mem=CALL_MEM
    }
    
    if (defined(TRUTH_VCF) && defined(TRUTH_VCF_INDEX)) {
    
        # To evaluate the VCF we need a template of the reference
        call buildReferenceTemplate {
            input:
                in_reference_file=reference_file
        }
        
        # Direct vcfeval comparison makes an archive with FP and FN VCFs
        call compareCalls {
            input:
                in_sample_vcf_file=bgzipMergedVCF.output_merged_vcf,
                in_sample_vcf_index_file=bgzipMergedVCF.output_merged_vcf_index,
                in_truth_vcf_file=select_first([TRUTH_VCF]),
                in_truth_vcf_index_file=select_first([TRUTH_VCF_INDEX]),
                in_template_archive=buildReferenceTemplate.output_template_archive,
                in_evaluation_regions_file=EVALUATION_REGIONS_BED,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
        
        # Hap.py comparison makes accuracy results stratified by SNPs and indels
        call compareCallsHappy {
            input:
                in_sample_vcf_file=bgzipMergedVCF.output_merged_vcf,
                in_sample_vcf_index_file=bgzipMergedVCF.output_merged_vcf_index,
                in_truth_vcf_file=select_first([TRUTH_VCF]),
                in_truth_vcf_index_file=select_first([TRUTH_VCF_INDEX]),
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_evaluation_regions_file=EVALUATION_REGIONS_BED,
                in_call_disk=CALL_DISK,
                in_call_mem=CALL_MEM
        }
    }
    
    output {
        File? output_vcfeval_evaluation_archive = compareCalls.output_evaluation_archive
        File? output_happy_evaluation_archive = compareCallsHappy.output_evaluation_archive
        File output_vcf = bgzipMergedVCF.output_merged_vcf
        File output_vcf_index = bgzipMergedVCF.output_merged_vcf_index
        Array[File] output_indelrealigned_bams = runAbraRealigner.indel_realigned_bam
        Array[File] output_indelrealigned_bam_indexes = runAbraRealigner.indel_realigned_bam_index
    }   
}

########################
### TASK DEFINITIONS ###
########################

task indexReference {
    input {
        File in_reference_file
        Int in_index_mem
        Int in_index_disk
    }

    command <<<
        set -eux -o pipefail
        
        ln -s ~{in_reference_file} ref.fa
        
        samtools faidx ref.fa
        
        # Save a reference copy by making the dict now
        java -jar /usr/picard/picard.jar CreateSequenceDictionary \
          R=ref.fa \
          O=ref.dict
    >>>
    output {
        File reference_index_file = "ref.fa.fai"
        File reference_dict_file = "ref.dict"
    }
    runtime {
        memory: in_index_mem + " GB"
        disks: "local-disk " + in_index_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task splitBAMbyPath {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        set -eux -o pipefail

        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai

        while read -r contig; do
            samtools view \
              -@ ~{in_map_cores} \
              -h -O BAM \
              input_bam_file.bam ${contig} \
              -o ~{in_sample_name}.${contig}.bam \
            && samtools index \
              ~{in_sample_name}.${contig}.bam
        done < "~{in_path_list_file}"
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai")
    }
    runtime {
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task runGATKRealignerTargetCreator {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_index_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_call_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command 
        # to exit with a non-zero status, or zero if all commands of the pipeline exit 
        set -o pipefail
        # cause a bash script to exit immediately when a command fails 
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately 
        set -u
        # echo each line of the script to stdout so we can see what is happening 
        set -o xtrace
        #to turn off echo do 'set +o xtrace' 

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s "~{in_reference_file}" reference.fa
        ln -f -s "~{in_reference_index_file}" reference.fa.fai
        # And the dict must be adjacent to both
        ln -f -s "~{in_reference_dict_file}" reference.dict

        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          --remove_program_records \
          -drf DuplicateRead \
          --disable_bam_indexing \
          -nt "32" \
          -R reference.fa \
          -L ${CONTIG_ID} \
          -I input_bam_file.bam \
          --out forIndelRealigner.intervals

        awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{in_sample_name}.${CONTIG_ID}.intervals.bed
    >>>
    output {
        File realigner_target_bed = glob("*.bed")[0]
    }
    runtime {
        time: 180
        memory: 20 + " GB"
        cpu: 16
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b"
    }
}
task widenRealignmentTargets {
    input {
        File in_target_bed_file
        File in_reference_index_file
        Int in_expansion_bases
        Int in_call_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        BASE_NAME=($(ls ~{in_target_bed_file} | rev | cut -f1 -d'/' | rev | sed s/.bed$//g))

        # Widen the BED regions, but don't escape the chromosomes
        bedtools slop -i "~{in_target_bed_file}" -g "~{in_reference_index_file}" -b "~{in_expansion_bases}" > "${BASE_NAME}.widened.bed"
    >>>
    output {
        File output_target_bed_file = glob("*.widened.bed")[0]
    }
    runtime {
        memory: 4 + " GB"
        cpu: 1
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    }
}
task runAbraRealigner {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_index_file
        File in_target_bed_file
        File in_reference_file
        File in_reference_index_file
        Int in_call_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai

        java -Xmx20G -jar /opt/abra2/abra2.jar \
          --targets ~{in_target_bed_file} \
          --in input_bam_file.bam \
          --out ~{in_sample_name}.${CONTIG_ID}.indel_realigned.bam \
          --ref reference.fa \
          --index \
          --threads 32
    >>>
    output {
        File indel_realigned_bam = glob("~{in_sample_name}.*.indel_realigned.bam")[0]
        File indel_realigned_bam_index = glob("~{in_sample_name}.*.indel_realigned.bai")[0]
    }
    runtime {
        time: 180
        memory: 20 + " GB"
        cpu: 16
        disks: "local-disk " + in_call_disk + " SSD"
        # This used to be docker: "dceoy/abra2:latest" but they moved the tag
        # and it stopped working. A known good version has been rehosted on
        # Quay in case Docker Hub deletes it.
        docker: "quay.io/adamnovak/dceoy-abra2@sha256:43d09d1c10220cfeab09e2763c2c5257884fa4457bcaa224f4e3796a28a24bba"
    }
}

task runDeepVariant {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_file_index
        File in_reference_file
        File in_reference_index_file
        File? in_model_meta_file
        File? in_model_index_file
        File? in_model_data_file
        Int in_min_mapq
        Int in_call_cores
        Int in_call_disk
        Int in_call_mem
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        ln -s ~{in_bam_file} input_bam_file.child.bam
        ln -s ~{in_bam_file_index} input_bam_file.child.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.indel_realigned.bam$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai
        
        MODEL_ARGS=()
        if [[ ! -z "~{in_model_meta_file}" ]] ; then
            # Model files must be adjacent and not at arbitrary paths
            ln -f -s "~{in_model_meta_file}" model.meta
            ln -f -s "~{in_model_index_file}" model.index
            ln -f -s "~{in_model_data_file}" model.data-00000-of-00001
            MODEL_ARGS=("--customized_model" "model")
        fi

        # When making examples, throw out any reads that are more likely to be
        # mismapped than not (MAPQ 3 or less)
        /opt/deepvariant/bin/run_deepvariant \
        --make_examples_extra_args 'min_mapping_quality=~{in_min_mapq}' \
        --model_type=WGS \
        --regions ${CONTIG_ID} \
        --ref=reference.fa \
        --reads=input_bam_file.child.bam \
        --output_vcf="~{in_sample_name}_deepvariant.vcf.gz" \
        --output_gvcf="~{in_sample_name}_deepvariant.g.vcf.gz" \
        --intermediate_results_dir=tmp_deepvariant \
        "${MODEL_ARGS[@]}" \
        --num_shards=16
    >>>
    output {
        File output_vcf_file = "~{in_sample_name}_deepvariant.vcf.gz"
        File output_gvcf_file = "~{in_sample_name}_deepvariant.g.vcf.gz"
    }
    runtime {
        memory: in_call_mem + " GB"
        cpu: 8
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
        nvidiaDriverVersion: "418.87.00"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "google/deepvariant:1.1.0-gpu"
    }
}

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
        Int in_call_disk
        Int in_call_mem
    }

    command {
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        for vcf_file in ${sep=" " in_clipped_vcf_chunk_files} ; do
            bcftools index "$vcf_file"
        done
        bcftools concat -a ${sep=" " in_clipped_vcf_chunk_files} | bcftools sort - > ${in_sample_name}_merged.vcf
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf"
    }
    runtime {
        time: 60
        memory: in_call_mem + " GB"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}

task bgzipMergedVCF {
    input {
        String in_sample_name
        File in_merged_vcf_file
        String in_vg_container
        Int in_call_disk
        Int in_call_mem
    }

    # TODO:
    #   If GVCF in in_merged_vcf_file then output_vcf_extension="gvcf" else output_vcf_extension="vcf"
    command {
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        bgzip -c ${in_merged_vcf_file} > ${in_sample_name}_merged.vcf.gz && \
        tabix -f -p vcf ${in_sample_name}_merged.vcf.gz
    }
    output {
        File output_merged_vcf = "${in_sample_name}_merged.vcf.gz"
        File output_merged_vcf_index = "${in_sample_name}_merged.vcf.gz.tbi"
    }
    runtime {
        time: 30
        memory: in_call_mem + " GB"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: in_vg_container
    }
}

task buildReferenceTemplate {
    input {
        File in_reference_file
    }
    command <<<
        set -eux -o pipefail
    
        rtg format -o template.sdf "~{in_reference_file}"
        tar -czf template.sdf.tar.gz template.sdf/
    >>>
    output {
        File output_template_archive = "template.sdf.tar.gz"
    }
    runtime {
        docker: "realtimegenomics/rtg-tools:3.12.1"
        memory: 4 + " GB"
        cpu: 1
        disks: "local-disk " + 10 + " SSD" 
    }
}

task compareCalls {
    input {
        File in_sample_vcf_file
        File in_sample_vcf_index_file
        File in_truth_vcf_file
        File in_truth_vcf_index_file
        File in_template_archive
        File? in_evaluation_regions_file
        Int in_call_disk
        Int in_call_mem
    }
    command <<<
        set -eux -o pipefail
    
        # Put sample and truth near their indexes
        ln -s "~{in_sample_vcf_file}" sample.vcf.gz
        ln -s "~{in_sample_vcf_index_file}" sample.vcf.gz.tbi
        ln -s "~{in_truth_vcf_file}" truth.vcf.gz
        ln -s "~{in_truth_vcf_index_file}" truth.vcf.gz.tbi
        
        # Set up template; we assume it drops a "template.sdf"
        tar -xf "~{in_template_archive}"
    
        rtg vcfeval \
            --baseline truth.vcf.gz \
            --calls sample.vcf.gz \
            ~{"--evaluation-regions=" + in_evaluation_regions_file} \
            --template template.sdf \
            --threads 32 \
            --output vcfeval_results
            
        tar -czf vcfeval_results.tar.gz vcfeval_results/
    >>>
    output {
        File output_evaluation_archive = "vcfeval_results.tar.gz"
    }
    runtime {
        docker: "realtimegenomics/rtg-tools:3.12.1"
        cpu: 32
        disks: "local-disk " + in_call_disk + " SSD"
        memory: in_call_mem + " GB"
    }
}

task compareCallsHappy {
    input {
        File in_sample_vcf_file
        File in_sample_vcf_index_file
        File in_truth_vcf_file
        File in_truth_vcf_index_file
        File in_reference_file
        File in_reference_index_file
        File? in_evaluation_regions_file
        Int in_call_disk
        Int in_call_mem
    }
    command <<<
        set -eux -o pipefail
    
        # Put sample and truth near their indexes
        ln -s "~{in_sample_vcf_file}" sample.vcf.gz
        ln -s "~{in_sample_vcf_index_file}" sample.vcf.gz.tbi
        ln -s "~{in_truth_vcf_file}" truth.vcf.gz
        ln -s "~{in_truth_vcf_index_file}" truth.vcf.gz.tbi
        
        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s "~{in_reference_file}" reference.fa
        ln -f -s "~{in_reference_index_file}" reference.fa.fai
        
        mkdir happy_results
   
        /opt/hap.py/bin/hap.py \
            truth.vcf.gz \
            sample.vcf.gz \
            ~{"-f " + in_evaluation_regions_file} \
            --reference reference.fa \
            --threads 32 \
            --engine=vcfeval \
            -o happy_results/eval
    
        tar -czf happy_results.tar.gz happy_results/
    >>>
    output {
        File output_evaluation_archive = "happy_results.tar.gz"
    }
    runtime {
        docker: "jmcdani20/hap.py:v0.3.12"
        cpu: 32
        disks: "local-disk " + in_call_disk + " SSD"
        memory: in_call_mem + " GB"
    }
}

