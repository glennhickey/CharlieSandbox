#!/bin/bash
cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals
module load python/3.7 singularity
#source /data/markellocj/test_vg_wdl_run/HG002_mapping_vcf_eval/toil_vg_test/toilvenv/bin/activate
source /data/markellocj/mapeval_tools/toilvg_construct_mapeval_venv/bin/activate
WORK_DIR="/data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG005_construct_haplotype_sim_graph_workdir"
#rm -fr ${WORK_DIR}/vg-construct-v3.2.2_grch38-hg005-outstore
#rm -fr ${WORK_DIR}/tmp-vg-construct-v3.2.2_grch38_hg005
#mkdir -p ${WORK_DIR}/vg-construct-v3.2.2_grch38-hg005-outstore
#mkdir -p ${WORK_DIR}/tmp-vg-construct-v3.2.2_grch38_hg005
export TMPDIR=${WORK_DIR}/tmp-vg-construct-v3.2.2_grch38_hg005
export XDG_RUNTIME_DIR=""
export TOIL_SLURM_ARGS="-t 20:00:00"
#toil clean ${WORK_DIR}/vg-construct-v3.2.2_grch38-hg005-jobstore

toil-vg construct \
--restart \
${WORK_DIR}/vg-construct-v3.2.2_grch38-hg005-jobstore \
${WORK_DIR}/vg-construct-v3.2.2_grch38-hg005-outstore \
--setEnv PATH=$PATH \
--batchSystem slurm \
--container Singularity \
--vg_docker quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
--realTimeLogging \
--logInfo \
--workDir ${WORK_DIR}/tmp-vg-construct-v3.2.2_grch38_hg005 \
--cleanWorkDir onSuccess \
--fasta /data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz \
--out_name hg005_sample_grch38 \
--logFile construct_baseline.v3.2.2_grch38_hg005.log \
--xg_index --gbwt_index \
--force_phasing True \
--pangenome \
--sample_graph HG005 \
--merge_graphs \
--keep_vcfs \
--fasta_regions \
--vcf /data/markellocj/benchmark_data/HG005_cohort/HG005_HG006_HG007_triophased.vcf.gz \
--vcf_phasing /data/markellocj/benchmark_data/HG005_cohort/HG005_HG006_HG007_triophased.vcf.gz \
--statePollingWait 120 \
--rescueJobsFrequency 120 \
--defaultMemory 100Gi

