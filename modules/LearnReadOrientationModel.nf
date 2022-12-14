#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2

container "${params.gatk4__container}"

process LearnReadOrientationModel {

        tag "LearnReadOrientationModel $bam_id"
        publishDir "$params.outdirB/", mode:'copy'



        input:
                path ('*')
		tuple val(bam_id) , file(bams)
		path base_path

        output:
                path ("${bam_id}-T_${bam_id}-N.read-orientation-model.tar.gz")


        shell:

        '''
        # Change this when running the complete ' mutect2' pipeline - followed by these filtering steps 

        ls !{base_path}/Somatic-ShortV/nextflow/make_PON_and_run_mutect2/Using_14SubIntervals_and_sarkMatching_gnomAD/results_mutect2/!{bam_id}*f1r2.*.tar.gz > !{bam_id}-T_!{bam_id}-N.unfiltered_f1r2.args


        gatk LearnReadOrientationModel --java-options "-Xmx58g" \
                --input !{bam_id}-T_!{bam_id}-N.unfiltered_f1r2.args \
                -O !{bam_id}-T_!{bam_id}-N.read-orientation-model.tar.gz



        '''


}

