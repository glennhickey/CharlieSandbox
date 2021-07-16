#!/bin/bash
module load bwa samtools/1.11 picard
COHORT_NAME="HG005"
SAMPLE_NAME="HG006"
LOW_COHORT_NAME=$(echo $COHORT_NAME | tr '[:upper:]' '[:lower:]')
LOW_SAMPLE_NAME=$(echo $SAMPLE_NAME | tr '[:upper:]' '[:lower:]')
cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem
time bwa mem \
-t 32 \
/data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna \
/data/markellocj/raw_read_data/${COHORT_NAME}_cohort/${SAMPLE_NAME}/${SAMPLE_NAME}.novaseq.pcr-free.30x.R1.fastq.gz \
/data/markellocj/raw_read_data/${COHORT_NAME}_cohort/${SAMPLE_NAME}/${SAMPLE_NAME}.novaseq.pcr-free.30x.R2.fastq.gz \
> bwamem_${LOW_SAMPLE_NAME}_wgs.sam

samtools sort -@ 32 -O BAM bwamem_${LOW_SAMPLE_NAME}_wgs.sam > bwamem_${LOW_SAMPLE_NAME}_wgs.positionsorted.bam && rm bwamem_${LOW_SAMPLE_NAME}_wgs.sam
samtools index -@ 32 bwamem_${LOW_SAMPLE_NAME}_wgs.positionsorted.bam

time java -Xmx100g -XX:ParallelGCThreads=32 -jar $PICARDJARPATH/picard.jar \
ReorderSam \
VALIDATION_STRINGENCY=SILENT \
INPUT=bwamem_${LOW_SAMPLE_NAME}_wgs.positionsorted.bam \
OUTPUT=bwamem_${LOW_SAMPLE_NAME}_wgs.positionsorted.reordered.bam \
SEQUENCE_DICTIONARY=/data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict

samtools index -@ 32 bwamem_${LOW_SAMPLE_NAME}_wgs.positionsorted.reordered.bam && rm bwamem_${LOW_SAMPLE_NAME}_wgs.positionsorted.bam bwamem_${LOW_SAMPLE_NAME}_wgs.positionsorted.bam.bai

mkdir bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    samtools view -@ 32 -b bwamem_${LOW_SAMPLE_NAME}_wgs.positionsorted.reordered.bam chr${CHR} > bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr/raw.chr${CHR}.bam
    samtools index -@ 32 bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr/raw.chr${CHR}.bam
done

rm indel_realignment_bwamem_grch38.${LOW_SAMPLE_NAME}.swarm
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    echo -e "module load singularity GATK/3.8-1 samtools/1.11 \
    && cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem/bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr/ \
    && samtools addreplacerg -@ 32 -O BAM -r ID:1 -r LB:lib1 -r SM:${SAMPLE_NAME} -r PL:illumina -r PU:unit1 raw.chr${CHR}.bam > raw.chr${CHR}.gatk_ready.bam \
    && samtools index -@ 32 raw.chr${CHR}.gatk_ready.bam \
    && java -jar \$GATK_JAR -T RealignerTargetCreator -R /data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna -L chr${CHR} -I raw.chr${CHR}.gatk_ready.bam -o ${CHR}.intervals \
    && awk -F '[:-]' 'BEGIN { OFS = \"\t\" } { if( \$3 == \"\") { print \$1, \$2-1, \$2 } else { print \$1, \$2-1, \$3}}' ${CHR}.intervals > ${CHR}.intervals.bed \
    && singularity run -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference docker://dceoy/abra2:latest --targets ${CHR}.intervals.bed --in raw.chr${CHR}.gatk_ready.bam --out raw.chr${CHR}.indel_realigned.bam --ref grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna --threads 16" >> indel_realignment_bwamem_grch38.${LOW_SAMPLE_NAME}.swarm
done

cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem
swarm_jid=$(swarm -f indel_realignment_bwamem_grch38.${LOW_SAMPLE_NAME}.swarm -t 16 -g 100 --partition=norm --time=12:00:00 --module singularity,GATK/3.8-1,samtools)

echo -e "cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem/bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr \
&& module load samtools/1.11 \
&& time samtools merge -@ 32 ${LOW_SAMPLE_NAME}.bwamem.abra_indel_realigned.merged.bam raw.chr{\?,\?\?}.indel_realigned.bam \
&& time samtools sort -@ 32 ${LOW_SAMPLE_NAME}.bwamem.abra_indel_realigned.merged.bam -O BAM > ${LOW_SAMPLE_NAME}.bwamem.abra_indel_realigned.positionsorted.bam && rm ${LOW_SAMPLE_NAME}.bwamem.abra_indel_realigned.merged.bam" >> merge_${LOW_SAMPLE_NAME}_indel_realigned_bams.sh

sbatch -c 32 --mem=100GB --time=32:00:00 --dependency=afterany:$swarm_jid merge_${LOW_SAMPLE_NAME}_indel_realigned_bams.sh

