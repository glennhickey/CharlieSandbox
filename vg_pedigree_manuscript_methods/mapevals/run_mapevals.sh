#!/bin/bash

# RUN BWAMEM
WORK_DIR="/data/Udpbinfo/usr/markellocj/vg_trio_methods/redo_mapevals/HG002_sim_reads"
#SIMULATED_READ_DIR="."
#SIMULATED_READ_DIR="high_conf_hg002_v4.2.1_regions"
#SIMULATED_READ_DIR="all_difficult_regions_hg002_v4.2.1_regions"
SIMULATED_READ_DIR="alllowmapandsegdupregions_hg002_v4.2.1_regions"

cd ${WORK_DIR}
module load samtools bwa singularity python
time bwa mem -t 32 -p /data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna ${WORK_DIR}/${SIMULATED_READ_DIR}/sim.paired.fq > mapped.bam

time samtools view -F 2048 -b mapped.bam > mapped.primary.bam
time samtools view -f 2048 -b mapped.bam > mapped.secondary.bam

singularity exec --env SIMULATED_READ_DIR=${SIMULATED_READ_DIR} -H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/vgteam/vg:v1.31.0 \
/bin/bash -c 'set -eo pipefail && \
SIMULATED_READ_DIR=${SIMULATED_READ_DIR} && \
time vg inject -x hg002_sample_grch38.xg mapped.primary.bam > mapped.primary.gam && \
time vg inject -x hg002_sample_grch38.xg mapped.secondary.bam > mapped.secondary.gam && \
time vg view -aj mapped.primary.gam | sed "s/\/1/_1/g" | sed "s/\/2/_2/g" | vg view -aGJ - | vg annotate -m -x hg002_sample_grch38.xg -a - | vg gamcompare -r 100 -s - ${SIMULATED_READ_DIR}/sim.gam 2> count | vg view -aj - > compared.primary.json && \
time vg view -aj mapped.secondary.gam | sed "s/\/1/_1/g" | sed "s/\/2/_2/g" | vg view -aGJ - | vg annotate -m -x hg002_sample_grch38.xg -a - | vg gamcompare -r 100 - ${SIMULATED_READ_DIR}/sim.gam | vg view -aj - > compared.secondary.json'

time python /data/markellocj/git/giraffe-sv-paper/scripts/mapping/combine_reads.py compared.primary.json compared.secondary.json bwamem.compared.json

singularity exec -H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/vgteam/vg:v1.31.0 \
/bin/bash -c 'set -eo pipefail && \
sed -i "/^$/d" bwamem.compared.json && \
printf "correct\tmq\tscore\taligner\n" > bwamem_roc_stats.tsv && \
READS="HG002" && \
PAIRING="paired" && \
GRAPH="hg002_sample_grch38" && \
ALGORITHM="bwamem" && \
CORRECT_COUNT="$(grep correctly_mapped bwamem.compared.json | wc -l)" && \
SCORE="$(sed -n "2p" count | sed "s/[^0-9\.]//g")" && \
MAPQ="$(grep mapping_quality\":\ 60 bwamem.compared.json | wc -l)" && \
MAPQ60="$(grep -v correctly_mapped bwamem.compared.json | grep mapping_quality\":\ 60 | wc -l)" && \
IDENTITY=$(jq ".identity" bwamem.compared.json | awk "{sum+=\$1} END {print sum/NR}") && \
echo ${GRAPH} ${READS} ${PAIRING} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${SCORE} && \
printf "graph\talgorithm\treads\tspeed\tpairing\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > report_bwamem.tsv && \
printf "${GRAPH}\t${ALGORITHM}\t${READS}\t${PAIRING}\t-\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> report_bwamem.tsv && \
jq -r "(if .correctly_mapped then 1 else 0 end|tostring) + \",\" + (.mapping_quality|tostring) + \",\" + (.score|tostring)" bwamem.compared.json | sed "s/\,/\t/g" | sed "s/$/\tbwamem_${GRAPH}${READS}${PAIRING}/" | sed "s/single//g ; s/paired/-pe/g ; s/null/0/g" >> bwamem_roc_stats.tsv'

# RUN GIRAFFE PRIMARY
singularity --env SIMULATED_READ_DIR=${SIMULATED_READ_DIR} -H ${PWD}:${HOME} --pwd ${HOME} \
-B /data/markellocj/graph_reference/snp1kg_ref_wgs/grch38_compact_decoys_primary_ref:${HOME}/grch38_compact_decoys_primary_ref docker://quay.io/vgteam/vg:v1.31.0 \
/bin/bash -c 'set -eo pipefail && \
PARAM_PRESET="default" && \
PAIRED="-i" && \
time vg giraffe -x grch38_compact_decoys_primary_ref/primarygrch38_no_alt_plus_hs38d1_v1.27.0_90_ga64b70c1f.xg -H grch38_compact_decoys_primary_ref/primarygrch38_no_alt_plus_hs38d1_v1.27.0_90_ga64b70c1f.gbwt -m grch38_compact_decoys_primary_ref/primarygrch38_no_alt_plus_hs38d1_v1.27.0_90_ga64b70c1f.min -d grch38_compact_decoys_primary_ref/primarygrch38_no_alt_plus_hs38d1_v1.27.0_90_ga64b70c1f.dist -G ${SIMULATED_READ_DIR}/sim.gam -b ${PARAM_PRESET} ${PAIRED} -t 22 > mapped.primary.gam && \
time vg gamcompare -r 100 -s <(vg annotate -m -x grch38_compact_decoys_primary_ref/primarygrch38_no_alt_plus_hs38d1_v1.27.0_90_ga64b70c1f.xg -a mapped.primary.gam) ${SIMULATED_READ_DIR}/sim.gam 2>count | vg view -aj - > giraffe_primary.compared.json && \
GRAPH="hg002_sample_grch38" && \
GBWT="sampled.64" && \
READS="HG002" && \
PARAM_PRESET="default" && \
PAIRING="paired" && \
CORRECT_COUNT="$(sed -n "1p" count | sed "s/[^0-9]//g")" && \
SCORE="$(sed -n "2p" count | sed "s/[^0-9\.]//g")" && \
MAPQ="$(grep mapping_quality\":\ 60 giraffe_primary.compared.json | wc -l)" && \
MAPQ60="$(grep -v correctly_mapped giraffe_primary.compared.json | grep mapping_quality\":\ 60 | wc -l)" && \
IDENTITY=$(jq ".identity" giraffe_primary.compared.json | awk "{sum+=\$1} END {print sum/NR}") && \
echo ${GRAPH} ${GBWT} ${READS} ${PARAM_PRESET}${PAIRING} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${IDENTITY} ${SCORE} && \
printf "graph\tgbwt\treads\tpairing\tspeed\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > giraffe_primary_report.tsv && \
printf "${GRAPH}\t${GBWT}\t${READS}\t${PARAM_PRESET}\t${PAIRING}\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> giraffe_primary_report.tsv && \
printf "correct\tmq\tscore\taligner\n" > giraffe_primary_roc_stats.tsv && \
jq -r "(if .correctly_mapped then 1 else 0 end|tostring) + \",\" + (.mapping_quality|tostring) + \",\" + (.score|tostring)" giraffe_primary.compared.json | sed "s/\,/\t/g" | sed "s/$/\tgiraffe_primary_${PARAM_PRESET}_${GRAPH}${GBWT}${READS}${PAIRING}/" >> giraffe_primary_roc_stats.tsv && \
sed -i "s/single//g ; s/paired/-pe/g ; s/null/0/g" giraffe_primary_roc_stats.tsv'

# RUN GIRAFFE 1000GP
singularity shell --env SIMULATED_READ_DIR=${SIMULATED_READ_DIR} -H ${PWD}:${HOME} --pwd ${HOME} \
-B /data/markellocj/graph_reference/snp1kg_ref_wgs/grch38_liftover_snp1kg_nosegdup/all_samples:${HOME}/all_samples docker://quay.io/vgteam/vg:v1.31.0 \
/bin/bash -c 'set -eo pipefail && \
PARAM_PRESET="default" && \
PAIRED="-i" && \
time vg giraffe -x all_samples/liftover_snp1kg_grch38_nosegdup.xg -H all_samples/liftover_snp1kg_grch38_nosegdup.gbwt -m all_samples/liftover_snp1kg_grch38_nosegdup.min -d all_samples/liftover_snp1kg_grch38_nosegdup.dist -G ${SIMULATED_READ_DIR}/sim.gam -b ${PARAM_PRESET} ${PAIRED} -t 22 > mapped.snp1kg.gam && \
time vg gamcompare -r 100 -s <(vg annotate -m -x all_samples/liftover_snp1kg_grch38_nosegdup.xg -a mapped.snp1kg.gam) ${SIMULATED_READ_DIR}/sim.gam 2>count | vg view -aj - > giraffe_snp1kg.compared.json && \
GRAPH="hg002_sample_grch38" && \
GBWT="sampled.64" && \
READS="HG002" && \
PARAM_PRESET="default" && \
PAIRING="paired" && \
CORRECT_COUNT="$(sed -n "1p" count | sed "s/[^0-9]//g")" && \
SCORE="$(sed -n "2p" count | sed "s/[^0-9\.]//g")" && \
MAPQ="$(grep mapping_quality\":\ 60 giraffe_snp1kg.compared.json | wc -l)" && \
MAPQ60="$(grep -v correctly_mapped giraffe_snp1kg.compared.json | grep mapping_quality\":\ 60 | wc -l)" && \
IDENTITY=$(jq ".identity" giraffe_snp1kg.compared.json | awk "{sum+=\$1} END {print sum/NR}") && \
echo ${GRAPH} ${GBWT} ${READS} ${PARAM_PRESET}${PAIRING} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${IDENTITY} ${SCORE} && \
printf "graph\tgbwt\treads\tpairing\tspeed\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > giraffe_snp1kg_report.tsv && \
printf "${GRAPH}\t${GBWT}\t${READS}\t${PARAM_PRESET}\t${PAIRING}\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> giraffe_snp1kg_report.tsv && \
printf "correct\tmq\tscore\taligner\n" > giraffe_snp1kg_roc_stats.tsv && \
jq -r "(if .correctly_mapped then 1 else 0 end|tostring) + \",\" + (.mapping_quality|tostring) + \",\" + (.score|tostring)" giraffe_snp1kg.compared.json | sed "s/\,/\t/g" | sed "s/$/\tgiraffe_snp1kg_${PARAM_PRESET}_${GRAPH}${GBWT}${READS}${PAIRING}/" >> giraffe_snp1kg_roc_stats.tsv && \
sed -i "s/single//g ; s/paired/-pe/g ; s/null/0/g" giraffe_snp1kg_roc_stats.tsv'

# RUN GIRAFFE PARENTAL
singularity shell --env SIMULATED_READ_DIR=${SIMULATED_READ_DIR} -H ${PWD}:${HOME} --pwd ${HOME} \
-B /data/markellocj/graph_reference/parental_ref_deeptrio_wgs_grch38_nosegdup_v1.31.0/hg002_cohort:${HOME}/hg002_cohort docker://quay.io/vgteam/vg:v1.31.0 \
/bin/bash -c 'set -eo pipefail && \
PARAM_PRESET="default" && \
PAIRED="-i" && \
time vg giraffe -x hg002_cohort/HG002.parental.graphs.xg -H hg002_cohort/HG002.parental.graphs.gbwt -m hg002_cohort/HG002.parental.graphs.min -d hg002_cohort/HG002.parental.graphs.dist -G ${SIMULATED_READ_DIR}/sim.gam -b ${PARAM_PRESET} ${PAIRED} -t 22 > mapped.parental.gam && \
time vg gamcompare -r 100 -s <(vg annotate -m -x hg002_cohort/HG002.parental.graphs.xg -a mapped.parental.gam) ${SIMULATED_READ_DIR}/sim.gam 2>count | vg view -aj - > giraffe_parental.compared.json && \
GRAPH="hg002_sample_grch38" && \
GBWT="sampled.64" && \
READS="HG002" && \
PARAM_PRESET="default" && \
PAIRING="paired" && \
CORRECT_COUNT="$(sed -n "1p" count | sed "s/[^0-9]//g")" && \
SCORE="$(sed -n "2p" count | sed "s/[^0-9\.]//g")" && \
MAPQ="$(grep mapping_quality\":\ 60 giraffe_parental.compared.json | wc -l)" && \
MAPQ60="$(grep -v correctly_mapped giraffe_parental.compared.json | grep mapping_quality\":\ 60 | wc -l)" && \
IDENTITY=$(jq ".identity" giraffe_parental.compared.json | awk "{sum+=\$1} END {print sum/NR}") && \
echo ${GRAPH} ${GBWT} ${READS} ${PARAM_PRESET}${PAIRING} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${IDENTITY} ${SCORE} && \
printf "graph\tgbwt\treads\tpairing\tspeed\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > giraffe_parental_report.tsv && \
printf "${GRAPH}\t${GBWT}\t${READS}\t${PARAM_PRESET}\t${PAIRING}\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> giraffe_parental_report.tsv && \
printf "correct\tmq\tscore\taligner\n" > giraffe_parental_roc_stats.tsv && \
jq -r "(if .correctly_mapped then 1 else 0 end|tostring) + \",\" + (.mapping_quality|tostring) + \",\" + (.score|tostring)" giraffe_parental.compared.json | sed "s/\,/\t/g" | sed "s/$/\tgiraffe_parental_${PARAM_PRESET}_${GRAPH}${GBWT}${READS}${PAIRING}/" >> giraffe_parental_roc_stats.tsv && \
sed -i "s/single//g ; s/paired/-pe/g ; s/null/0/g" giraffe_parental_roc_stats.tsv'

