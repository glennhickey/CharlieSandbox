#!/bin/bash
module load singularity
cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG005_sim_reads

## SIMULATE BASELINE READS
#singularity shell -H ${PWD}:${HOME} --pwd ${HOME} \
#-B /data/markellocj/raw_read_data/HG005_cohort_precision_fda_reads/HG005:${HOME}/HG005 \
#docker://quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975
#
#NREADS=1000000
#FASTQ=HG005/HG005_merged_interleaved-shuffled-1m.fastq.gz
#
#echo simulate reads
#vg sim -r -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x hg005_sample_grch38.xg -g hg005_sample_grch38.gbwt --sample-name HG005 --ploidy-regex "hs38d1:0,chrNC_007605:0,chrX:1,chrY:1,chrY_.*:1,chrEBV:0,.*:2" -F $FASTQ > sim.raw.gam

declare -a REGION_LIST=( "high_conf_hg005_v4.2.1_regions" "all_difficult_regions_hg005_v4.2.1_regions" "alllowmapandsegdupregions_hg005_v4.2.1_regions" "mhc_hg005_v4.2.1_regions" "high_conf_NO1000GP_hg005_v4.2.1_regions" )
declare -a BED_FILE_LIST=( "HG005_GRCh38_1_22_draft_2_v4.2.1_benchmark.bed" "HG005_GRCh38_v4.2.1.all_difficult_regions.bed" "HG005_GRCh38_v4.2.1.alllowmapandsegdupregions.bed" "HG005_GRCh38_v4.2.1.MHC.bed" "HG005_GRCh38_CHROM1-22_v4.2.1.highconf.NO_SNP1KG.bed" )

for index in "${!REGION_LIST[@]}"; do
    cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG005_sim_reads
    REGION="${REGION_LIST[index]}"
    BED_FILE="${BED_FILE_LIST[index]}"
    mkdir -p ${REGION}
    cp sim.raw.gam hg005_sample_grch38.xg ${REGION}/
    cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG005_sim_reads/${REGION}
    singularity exec --env REGION=${REGION},BED_FILE=${BED_FILE} -H ${PWD}:${HOME} --pwd ${HOME} \
    -B /data/markellocj/benchmark_data/HG005_cohort:${HOME}/HG005_cohort \
    -B /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG005_construct_haplotype_sim_graph_workdir/vg-construct-v3.2.2_grch38-hg005-outstore:${HOME}/vg-construct-v3.2.2_grch38-hg005-outstore \
    docker://quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
    /bin/bash -xc 'set -eo pipefail && \
    echo annotate with paths && \
    time vg annotate \
    -p -x vg-construct-v3.2.2_grch38-hg005-outstore/hg005_sample_grch38.vg -a sim.raw.gam \
    --bed-name HG005_cohort/${BED_FILE} \
    --threads 32 > sim.${REGION}.gam && \
    time vg filter -i -U -F "" sim.${REGION}.gam > sim.gam && \
    rm -f sim.${REGION}.gam sim.fq.gz sim.fq && \
    echo convert to fastq && \
    time vg view -X -a sim.gam | gzip > sim.fq.gz && \
    gunzip sim.fq.gz && \
    sed "s/_1\$//g" sim.fq | sed "s/_2\$//g" > sim.paired.fq && \
    echo format true position information && \
    time vg view -a sim.gam | jq -c -r "[.name] + if (.annotation.features | length) > 0 then [.annotation.features | join(\",\")] else [\".\"] end + if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + [.score] + if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv" > true.pos'
done




