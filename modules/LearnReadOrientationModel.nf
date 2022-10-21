#!/usr/bin/env nextflow

refdir='/scratch/wz54/gs5517/sarek_testing/Reference/v0'
base_path='/scratch/wz54/npd561/PIPE-2629_thyroid_carcinoma/nextflow_pipelines/nfcore/sarek'

process LearnReadOrientationModel {

        tag "LearnReadOrientationModel $bam_id"
        publishDir "$params.outdirB/", mode:'copy'



        input:
                path ('*')
		tuple val(bam_id) , file(bams)


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

