#!/bin/bash
cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/hg002_cohort/bwamem/;
bcftools concat -O v HG002.bwamem.deeptrio.chr{?,??}.vcf.gz > merged.deeptrio.vcf && bcftools sort merged.deeptrio.vcf -O v > merged.deeptrio.sorted.vcf && rm merged.deeptrio.vcf && bgzip merged.deeptrio.sorted.vcf && tabix -p vcf merged.deeptrio.sorted.vcf.gz
