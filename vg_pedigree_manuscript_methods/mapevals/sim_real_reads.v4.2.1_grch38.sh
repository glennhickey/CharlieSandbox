#!/bin/bash
module load python/3.7 singularity
cd /data/Udpbinfo/usr/markellocj/sim_map_eval_150bp
source /data/markellocj/test_toil_vg_run/toil_vg_tools/toilvg_venv/bin/activate
WORK_DIR="/data/Udpbinfo/usr/markellocj/sim_map_eval_150bp/simulate_reads_workdir"
rm -fr ${WORK_DIR}/vg-sim-v4.2.1_grch38-outstore
rm -fr ${WORK_DIR}/tmp-vg-sim-v4.2.1_grch38
mkdir -p ${WORK_DIR}/vg-sim-v4.2.1_grch38-outstore
mkdir -p ${WORK_DIR}/tmp-vg-sim-v4.2.1_grch38
cd ${WORK_DIR}
export TMPDIR=${WORK_DIR}/tmp-vg-sim-v4.2.1_grch38
toil clean ${WORK_DIR}/vg-sim-v4.2.1_grch38-jobstore
toil-vg sim \
${WORK_DIR}/vg-sim-v4.2.1_grch38-jobstore \
/data/Udpbinfo/usr/markellocj/sim_map_eval_150bp/construct_haplotype_graph_workdir/vg-construct-v4.2.1_grch38-outstore/baseline_HG002_haplo_thread_0.xg \
/data/Udpbinfo/usr/markellocj/sim_map_eval_150bp/construct_haplotype_graph_workdir/vg-construct-v4.2.1_grch38-outstore/baseline_HG002_haplo_thread_1.xg \
65000000 \
${WORK_DIR}/vg-sim-v4.2.1_grch38-outstore \
--batchSystem singleMachine \
--container Singularity \
--logFile sim_baseline_HG002_haplo_sim_65000000_ekg-s1.hg002_150bp_realreads.v4.2.1_grch38.log \
--out_name baseline_HG002_haplo_sim_65000000_ekg-s1_hg002_150bp_realreads_v4.2.1_grch38 \
--sim_chunks 1 \
--gam --maxCores 32 \
--sim_opts ' -I -p 500 -v 50' --seed 1 \
--fastq /data/markellocj/raw_read_data/HG002_cohort_precision_fda_reads/HG002/HG002_merged_interleaved.fastq.gz \
--fastq_out \
--whole_genome_config \
--annotate_xg /data/Udpbinfo/usr/markellocj/sim_map_eval_150bp/construct_haplotype_graph_workdir/vg-construct-v4.2.1_grch38-outstore/baseline_HG002_haplo.xg
