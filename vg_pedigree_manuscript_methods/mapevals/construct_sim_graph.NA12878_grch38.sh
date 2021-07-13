#!/bin/bash
cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals
module load python/3.7 singularity
source /data/markellocj/test_toil_vg_run/toil_vg_tools/toilvg_venv/bin/activate
WORK_DIR="/data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/NA12878_construct_haplotype_sim_graph_workdir"
rm -fr ${WORK_DIR}/vg-construct_grch38-outstore
rm -fr ${WORK_DIR}/tmp-vg-construct_grch38
mkdir -p ${WORK_DIR}/vg-construct_grch38-outstore
mkdir -p ${WORK_DIR}/tmp-vg-construct_grch38
toil clean ${WORK_DIR}/vg-construct_grch38-jobstore

VCFS=()
## NOT including chrY because there are no genotype calls for NA12878 in the 1000GP VCF
for CHROM in {1..22} X ; do
  VCFS+=("/data/markellocj/graph_reference/construction_source_files/1kg_phased_genotypes_grch38/liftover_from_grch37_no_segdubs_gt10kb/ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.renamed.vcf.gz")
done

toil-vg construct \
${WORK_DIR}/vg-construct_grch38-jobstore \
${WORK_DIR}/vg-construct_grch38-outstore \
--batchSystem singleMachine \
--container Singularity \
--vg_docker quay.io/vgteam/vg:v1.31.0 \
--realTimeLogging \
--logInfo \
--workDir ${WORK_DIR}/tmp-vg-construct_grch38 \
--cleanWorkDir onSuccess \
--fasta /data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz \
--out_name hg002_sample_grch38 \
--logFile construct_baseline_grch38.log \
--xg_index --gbwt_index \
--force_phasing True \
--sample_graph NA12878 \
--merge_graphs \
--keep_vcfs \
--fasta_regions \
--vcf "${VCFS[@]}" \
--vcf_phasing "${VCFS[@]}" \
--statePollingWait 120 \
--rescueJobsFrequency 120

