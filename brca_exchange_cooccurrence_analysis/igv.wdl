version 1.0

## Copyright Broad Institute, 2020
##
## This WDL defines the workflow for generating a tarball containing IGV screenshots.
## Each screenshot contains tracks for proband WES, father, mother, and proband WGS.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

###########################################################################
# WORKFLOW DEFINITION
###########################################################################

workflow generate_igv_screenshots {
    input {

        File batch_script # make_igv_batchfile.py

        File input_table

        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict

        String output_folder_prefix 

    }

    ##########################################################################
    ## for each variant in input_table, 
    ##     1. run GATK PrintReads() to generate mini-crams for pb_WES, fa, mo, pb_WGS
    ##     2. run IGV in batch mode to produce single screenshot
    ##########################################################################

    ## Read in table and coerce as Array[String] corresponding to a sample
    ## 4-column format
    ##      sample_id, var_id, region, sample cram google bucket path
    call read_table {
        input:
            table = input_table
    }

    ## i index (rows) correspond to individual samples
    scatter (i in range(length(read_table.out))) {

        String sample_id = read_table.out[i][0]
        String var_ids = read_table.out[i][1]
        String regions = read_table.out[i][2]

        # note: keep these as strings to avoid localizing large cram files (GATK will stream from bucket)
        String pb_cram = read_table.out[i][3]

        call generate_mini_crams {
            input:
                var_id = var_ids,
                region = regions,
                pb_cram = pb_cram,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_fasta_dict = ref_fasta_dict
        }

        call run_igv {
            input:
                sample_id = sample_id,
                var_id = var_ids,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                script = batch_script,
                pb_bam = generate_mini_crams.mini_pb_bam,
                pb_bam_index = generate_mini_crams.mini_pb_bam_index
        }

    }

    call gather_shards {
        input:
            screens = run_igv.screenshots,
            outprefix = output_folder_prefix

    }

    output {
        File screenshots_folder = gather_shards.screenshots_tarball
    }

}


###########################################################################
# TASK DEFINITIONS
###########################################################################


## Note: requires specific 4-column format:
## 4-column format
##      sample_id, var_id, region, sample cram google bucket path
task read_table {
    input {
        File table
    }

    command { 
        echo "reading table" 
    }

    runtime {
        docker: "ubuntu:latest"
        preemptible: 3
        maxRetries: 3
    }

    output{
        Array[Array[String]] out = read_tsv(table)
    }

}


## Uses GATK PrintReads to generate mini cram files for pb WES, fa, mo, pb WGS
task generate_mini_crams {
    input {

        String var_id

        String region
        
        String pb_cram

        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict

        Int disk_size = 50

        String project_id  # test if this google project ID works

    }
    #String region_list = `echo ~{region} | tr , \\n`


    command <<< 

        echo ~{region} | tr , \\n > tmp.region.list

        ## how to set project for --gcs-project-for-requester-pays option??
        gatk PrintReads -I ~{pb_cram} -L tmp.region.list -O "pb_wes.bam" -R ~{ref_fasta} --gcs-project-for-requester-pays ~{project_id}
    >>>

    runtime {
        #docker: "broadinstitute/gatk:4.1.8.1"
        docker: "broadinstitute/gatk:latest"
        memory: "8G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        maxRetries: 3
    }

    output{
        File region_list = "tmp.region.list"
        File mini_pb_bam = "pb_wes.bam"
        File mini_pb_bam_index = "pb_wes.bai"
    }

}


## Run IGV to produce screenshots for 500-bp surrounding each variant
## 4 tracks: sample WES, father, mother, sample WGS
task run_igv {
    input {
        String sample_id
        String var_id
        File ref_fasta
        File ref_fasta_index
        
        File script
        
        File pb_bam
        File pb_bam_index 
        
        Int disk_size = 100
    }

    command { 
        
        ## generate IGV batch file + prints out screenshot filename and stores in bash variable $SCREENSHOT
        python ~{script} -s ~{sample_id} -v ~{var_id} -r ~{ref_fasta} -b ~{pb_bam} -o batch.txt
        

        ## RUN IGV IN BATCH MODE
        xvfb-run --server-args="-screen 0, 1920x540x24" bash /IGV_2.4.14/igv.sh -b batch.txt

        ## create output directory for screenshots
        mkdir screenshots

        ## move screens to new folder
        cp *.png screenshots

        ## create tarball of screenshots directory
        tar -czf screenshots.tar.gz screenshots
    }

    runtime {
        #docker: "talkowski/igv_gatk:latest"
        docker: "talkowski/igv_gatk"
        memory: "10G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        maxRetries: 3
    }

    output {
        File batchfile = "batch.txt"
        File screenshots = "screenshots.tar.gz"
    }
}


## Combines screenshot tarball outputs from each shard into single tarball
task gather_shards {
    input {
        Array[File] screens
        String outprefix
    }

    command <<<

        echo "creating directory ~{outprefix}/ ..."

        mkdir ~{outprefix}

        echo "uncompressing tarballs and copying screenshots into ~{outprefix}/ ..."
        while read file; do
        echo "uncompressing ${file} ..."
        tar -zxf ${file}

        echo "moving screenshots"
        mv screenshots/* ~{outprefix}/

        echo "cleaning up directories"
        rm -rf ${file} screenshot/

        done < ~{write_lines(screens)};

        echo "creating ~{outprefix}.tar.gz ..."
        tar -czf ~{outprefix}.tar.gz ~{outprefix}/

        echo "done!"
    >>>

    runtime {
        docker: "ubuntu:latest"
        preemptible: 3
        maxRetries: 3
    }

    output {
        File screenshots_tarball = "~{outprefix}.tar.gz"
    }   

}
