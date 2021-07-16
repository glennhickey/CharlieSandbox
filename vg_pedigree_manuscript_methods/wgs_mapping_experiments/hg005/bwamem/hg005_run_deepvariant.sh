#!/bin/bash
module load bwa samtools/1.11 picard bcftools singularity
COHORT_NAME="HG005"
SAMPLE_NAME="HG005"
LOW_COHORT_NAME=$(echo $COHORT_NAME | tr '[:upper:]' '[:lower:]')
LOW_SAMPLE_NAME=$(echo $SAMPLE_NAME | tr '[:upper:]' '[:lower:]')
cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem/bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    samtools index -@ 32 raw.chr${CHR}.indel_realigned.bam
done

rm deepvariant_calling.${LOW_SAMPLE_NAME}.default.minmapq1.indel_realigned.swarm
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    echo -e "cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem/bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr; singularity run --nv -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference docker://google/deepvariant:1.1.0-gpu /opt/deepvariant/bin/run_deepvariant --make_examples_extra_args 'min_mapping_quality=1' --model_type=WGS --regions chr${CHR} --ref=grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna --reads=raw.chr${CHR}.indel_realigned.bam --output_vcf=chr${CHR}.indel_realigned.deepvariant.vcf.gz --output_gvcf=chr${CHR}.indel_realigned.deepvariant.g.vcf.gz --intermediate_results_dir=tmp_deepvariant.indel_realigned.chr${CHR} --num_shards=16" >> deepvariant_calling.${LOW_SAMPLE_NAME}.default.minmapq1.indel_realigned.swarm
done

swarm_jid=$(swarm -f deepvariant_calling.${LOW_SAMPLE_NAME}.default.minmapq1.indel_realigned.swarm -t 16 -g 20 --partition=gpu --gres=gpu:k80:1 --time=3:00:00 --module singularity)

echo -e "#!/bin/bash
cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem/bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr \
bcftools concat -O v chr{?,??}.indel_realigned.deepvariant.vcf.gz > merged.deepvariant.indel_realigned.vcf \
&& bcftools sort merged.deepvariant.indel_realigned.vcf -O v > merged.deepvariant.indel_realigned.sorted.vcf && rm merged.deepvariant.indel_realigned.vcf \
&& bgzip merged.deepvariant.indel_realigned.DEBUG.sorted.vcf \
&& tabix -p vcf merged.deepvariant.indel_realigned.DEBUG.sorted.vcf.gz" >> merge_${LOW_SAMPLE_NAME}_deepvariant_vcfs.sh

sbatch -c 8 --mem=50GB --time=4:00:00 --dependency=afterany:$swarm_jid merge_${LOW_SAMPLE_NAME}_deepvariant_vcfs.sh

