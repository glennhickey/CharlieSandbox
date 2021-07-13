#!/bin/bash
module load singularity
cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG002_sim_reads
singularity shell -H ${PWD}:${HOME} --pwd ${HOME} \
-B /data/markellocj/raw_read_data/HG002_cohort_precision_fda_reads/HG002:${HOME}/HG002 \
docker://quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975

NREADS=1000000
FASTQ=HG002/HG002_merged_interleaved-shuffled-1m.fastq.gz

echo simulate reads
vg sim -r -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x hg002_sample_grch38.xg -g hg002_sample_grch38.gbwt --sample-name HG002 --ploidy-regex "hs38d1:0,chrNC_007605:0,chrX:1,chrY:1,chrY_.*:1,chrEBV:0,.*:2" -F $FASTQ > sim.raw.gam

echo annotate with paths
vg annotate -p -x graph.vg -a sim.raw.gam > sim.gam

echo convert to fastq
vg view -X -a sim.gam | gzip > sim.fq.gz

echo format true position information
vg view -a sim.gam | jq -c -r '[.name] + if (.annotation.features | length) > 0 then [.annotation.features | join(",")] else ["."] end + if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + [.score] + if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv' > true.pos


