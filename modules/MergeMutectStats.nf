#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2


container "${params.gatk4__container}"

process MergeMutectStats {

        tag "MergeMutectStats"
        publishDir "$params.outdirB/", mode:'copy'


        input:
		tuple val(bam_id) , file(bams)
                path ('*') 
		path base_path


        output:
                path ("${bam_id}-T_${bam_id}-N.unfiltered_stats.args")
                path ("${bam_id}-T_${bam_id}-N.unfiltered.stats")


        shell:

        '''

        # Change this when running the complete ' mutect2' pipeline - followed by these filtering steps         
        #ls !{base_path}/Somatic-ShortV/nextflow/make_PON_and_run_mutect2/Using_14SubIntervals_and_sarkMatching_gnomAD/results_mutect2/!{bam_id}*.stats   >!{bam_id}-T_!{bam_id}-N.unfiltered_stats.args

        ls !{params.outdirA}/!{bam_id}*.stats   >!{bam_id}-T_!{bam_id}-N.unfiltered_stats.args



        gatk MergeMutectStats \
                --stats !{bam_id}-T_!{bam_id}-N.unfiltered_stats.args \
                -O !{bam_id}-T_!{bam_id}-N.unfiltered.stats



        '''


        }


