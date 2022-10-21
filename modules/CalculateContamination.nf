#!/usr/bin/env nextflow

refdir='/scratch/wz54/gs5517/sarek_testing/Reference/v0'
base_path='/scratch/wz54/npd561/PIPE-2629_thyroid_carcinoma/nextflow_pipelines/nfcore/sarek'

process CalculateContamination {

        tag "CalculateContamination $bam_id"
        publishDir "$params.outdirB/", mode:'copy'




        input:
		tuple val(bam_id) , file(bams)
                path pileupsTable_T
                path pileupsTable_N

                path patch


        output:
                path ("${bam_id}-T_${bam_id}-N_contamination.table")
                path ("${bam_id}-T_segments.table")

        shell:
        '''

        gatk  --java-options "-Xmx16g" \
                CalculateContamination \
                -I !{bam_id}-T.pileups.table \
                -tumor-segmentation !{bam_id}-T_segments.table \
                -matched !{bam_id}-N.pileups.table \
                -O !{bam_id}-T_!{bam_id}-N_contamination.table
        '''

        }


