#!/usr/bin/env nextflow

refdir='/scratch/wz54/gs5517/sarek_testing/Reference/v0'
base_path='/scratch/wz54/npd561/PIPE-2629_thyroid_carcinoma/nextflow_pipelines/nfcore/sarek'


process GatherVcfs_step {

        tag "GatherVcfs_step $bam_id"
        publishDir "$params.outdirA/Mutect2", mode:'copy'


        input:
                path ('*')
                tuple val(bam_id) , file(bams)

        output:
                path ("${bam_id}-T_${bam_id}-N.unfiltered.vcf.gz")
                path ("${bam_id}-T_${bam_id}-N.unfiltered.vcf.gz.tbi")



        script:

        """
        #ls $base_path/Somatic-ShortV/nextflow/make_PON_and_run_mutect2/final_scripts_runs_DSL2/results_mutect2/${bam_id}-T_${bam_id}-N.unfiltered.*.vcf.gz   >${bam_id}_gathered_vcfs_across_subintervals.list

        ls ${params.outdirA}/${bam_id}-T_${bam_id}-N.unfiltered.*.vcf.gz   >${bam_id}_gathered_vcfs_across_subintervals.list

        # GatherVcfs requires intervals in order, so add chrM using MergeVcfs
        gatk GatherVcfs \
                -I  ${bam_id}_gathered_vcfs_across_subintervals.list \
                -O  ${bam_id}-T_${bam_id}-N.unfiltered_unsorted.vcf.gz

        #Sort
        gatk SortVcf \
                -I ${bam_id}-T_${bam_id}-N.unfiltered_unsorted.vcf.gz \
                -O ${bam_id}-T_${bam_id}-N.unfiltered.vcf.gz
        """

}

