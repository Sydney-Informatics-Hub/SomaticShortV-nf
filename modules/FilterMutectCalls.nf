#!/usr/bin/env nextflow

refdir='/scratch/wz54/gs5517/sarek_testing/Reference/v0'
base_path='/scratch/wz54/npd561/PIPE-2629_thyroid_carcinoma/nextflow_pipelines/nfcore/sarek'

process FilterMutectCalls {

        tag "FilterMutectCalls $bam_id"
        publishDir "$params.outdirB/filtered/", mode:'copy'


        input :
		tuple val(bam_id) , file(bams)
                path pair_contaminationTable

                path outDirA
                path outdirB

        output :
                path ("${bam_id}-T_${bam_id}-N.filtered.vcf.gz")

        shell:
        '''
        gatk --java-options "-Xmx16g -Xms16g" \
                FilterMutectCalls \
                --reference !{refdir}/Homo_sapiens_assembly38.fasta \
                -V !{params.outdirA}/Mutect2/!{bam_id}-T_!{bam_id}-N.unfiltered.vcf.gz \
                --stats !{params.outdirB}/!{bam_id}-T_!{bam_id}-N.unfiltered.stats \
                --tumor-segmentation !{params.outdirB}/!{bam_id}-T_segments.table \
                --contamination-table  !{params.outdirB}/!{bam_id}-T_!{bam_id}-N_contamination.table \
                --ob-priors !{params.outdirB}/!{bam_id}-T_!{bam_id}-N.read-orientation-model.tar.gz \
                -O !{bam_id}-T_!{bam_id}-N.filtered.vcf.gz



        '''


        }

