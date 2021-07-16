#!/bin/bash
module load singularity python/3.7
source /data/Udpbinfo/usr/markellocj/test_vg_pedigree_giraffe_grch38/tiny_unittest_workdir/toilvg_test_venv/bin/activate
GRAPH_REF_DIR="/data/markellocj/graph_reference/snp1kg_ref_wgs/grch38_liftover_snp1kg_nosegdup/all_samples"
INPUT_DIR="/data/Udpbinfo/Scratch/markellocj/toil_vg_workflow_inputs/toil_vg_inputs/grch38_inputs"
WORKDIR="/data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/hg005_cohort"
OUTSTORE="${WORKDIR}/HG005_pedigree_giraffe_deeptrio_grch38_liftover_non_segdup_snp1kg_wgs_vg_1.31.0_outstore"
JOBSTORE="${WORKDIR}/HG005_pedigree_giraffe_deeptrio_grch38_liftover_non_segdup_snp1kg_wgs_vg_1.31.0_jobstore"
LOGFILE="${WORKDIR}/HG005_pedigree_giraffe_deeptrio_grch38_liftover_non_segdup_snp1kg_wgs_vg_1.31.0.log"
TMPDIR="${WORKDIR}/tmp_dev_grch38_liftover_vg_1.31.0"
export TOIL_SLURM_ARGS='-t 20:00:00'
export SINGULARITY_CACHEDIR=/data/markellocj/singularity_cache
#rm -fr ${LOGFILE} ${OUTSTORE} ${TMPDIR}
mkdir -p ${OUTSTORE} ${TMPDIR}
cd ${WORKDIR}

#toil clean ${JOBSTORE}

time toil-vg pedigree \
--restart \
--retryCount 0 \
--rotatingLogging \
--setEnv PATH=$PATH \
--disableProgress \
--realTimeLogging \
--batchSystem slurm \
--statePollingWait 120 \
--rescueJobsFrequency 120 \
--container Singularity \
--logInfo \
--logFile ${LOGFILE} \
--workDir ${TMPDIR} \
--cleanWorkDir onSuccess \
--whole_genome_config \
--vg_docker 'quay.io/vgteam/vg:v1.31.0' \
${JOBSTORE} \
${OUTSTORE} \
HG005 \
HG007 \
HG006 \
--fastq_proband /data/markellocj/raw_read_data/HG005_cohort/HG005/HG005.novaseq.pcr-free.30x.R1.fastq.gz /data/markellocj/raw_read_data/HG005_cohort/HG005/HG005.novaseq.pcr-free.30x.R2.fastq.gz \
--fastq_maternal /data/markellocj/raw_read_data/HG005_cohort/HG007/HG007.novaseq.pcr-free.30x.R1.fastq.gz /data/markellocj/raw_read_data/HG005_cohort/HG007/HG007.novaseq.pcr-free.30x.R2.fastq.gz \
--fastq_paternal /data/markellocj/raw_read_data/HG005_cohort/HG006/HG006.novaseq.pcr-free.30x.R1.fastq.gz /data/markellocj/raw_read_data/HG005_cohort/HG006/HG006.novaseq.pcr-free.30x.R2.fastq.gz \
--ref_fasta ${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna \
--ref_fasta_index ${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai \
--ref_fasta_dict ${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict \
--caller deepvariant \
--mapper giraffe \
--use_haplotypes \
--xg_index ${GRAPH_REF_DIR}/liftover_snp1kg_grch38_nosegdup.xg \
--gbwt_index ${GRAPH_REF_DIR}/liftover_snp1kg_grch38_nosegdup.gbwt \
--graph_gbwt_index ${GRAPH_REF_DIR}/liftover_snp1kg_grch38_nosegdup.gg \
--minimizer_index ${GRAPH_REF_DIR}/liftover_snp1kg_grch38_nosegdup.min \
--distance_index ${GRAPH_REF_DIR}/liftover_snp1kg_grch38_nosegdup.dist \
--id_ranges ${INPUT_DIR}/path_list_whole_genome.txt \
--path_list ${INPUT_DIR}/path_list_whole_genome.txt \
--ped_file /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/hg005_cohort/HG005.ped \
--eagle_data ${INPUT_DIR}/eagle_data_grch38.tar.gz \
--bam_output \
--use_decoys 2> stdout.debug.txt


