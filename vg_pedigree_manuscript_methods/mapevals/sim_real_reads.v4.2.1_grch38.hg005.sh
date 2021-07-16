#!/bin/bash

module load bbtools seqtk
cd /data/markellocj/raw_read_data/HG005_cohort/HG005/
zcat HG005.novaseq.pcr-free.30x.R1.fastq.gz | seqtk sample -s100 - 1000000 | gzip - > HG005.novaseq.pcr-free.30x.R1-1m.fq.gz
zcat HG005.novaseq.pcr-free.30x.R2.fastq.gz | seqtk sample -s100 - 1000000 | gzip - > HG005.novaseq.pcr-free.30x.R2-1m.fq.gz
paste <(zcat HG005.novaseq.pcr-free.30x.R1-1m.fq.gz) <(zcat HG005.novaseq.pcr-free.30x.R2-1m.fq.gz) | paste - - - - | shuf | awk -F'\t' '{OFS="\n"; print $1,$3,$5,$7 > "HG005.novaseq.pcr-free.30x.R1-shuffled-1m.fq"; print $2,$4,$6,$8 > "HG005.novaseq.pcr-free.30x.R2-shuffled-1m.fq"}'
gzip HG005.novaseq.pcr-free.30x.R1-shuffled-1m.fq
gzip HG005.novaseq.pcr-free.30x.R2-shuffled-1m.fq
bbtools reformat t=32 in=HG005.novaseq.pcr-free.30x.R1-shuffled-1m.fq.gz in2=HG005.novaseq.pcr-free.30x.R2-shuffled-1m.fq.gz out=HG005_merged_interleaved-shuffled-1m.fastq.gz

mkdir -p /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG005_sim_reads
cd /data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG005_sim_reads
singularity shell -H ${PWD}:${HOME} --pwd ${HOME} \
-B /data/markellocj/raw_read_data/HG005_cohort/HG005:${HOME}/HG005 \
docker://quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975

NREADS=1000000
FASTQ=HG005/HG005_merged_interleaved-shuffled-1m.fastq.gz

echo simulate reads
time vg sim -r -I -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x hg005_sample_grch38.xg -g hg005_sample_grch38.gbwt --sample-name HG005 --ploidy-regex "hs38d1:0,chrNC_007605:0,chrX:1,chrY:1,chrY_.*:1,chrEBV:0,.*:2" -F $FASTQ > sim.raw.gam

echo annotate with paths
time vg annotate -p -x hg005_sample_grch38.vg -a sim.raw.gam > sim.gam

echo convert to fastq
time vg view -X -a sim.gam | gzip > sim.fq.gz

gunzip sim.fq.gz
sed 's/_1$//g' sim.fq | sed 's/_2$//g' > sim.paired.fq

echo format true position information
time vg view -a sim.gam | jq -c -r '[.name] + if (.annotation.features | length) > 0 then [.annotation.features | join(",")] else ["."] end + if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + [.score] + if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv' > true.pos


