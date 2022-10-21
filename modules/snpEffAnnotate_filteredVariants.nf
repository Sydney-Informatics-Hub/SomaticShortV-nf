#!/usr/bin/env nextflow

refdir='/scratch/wz54/gs5517/sarek_testing/Reference/v0'
base_path='/scratch/wz54/npd561/PIPE-2629_thyroid_carcinoma/nextflow_pipelines/nfcore/sarek'

params.outdirB="$base_path/Somatic-ShortV/nextflow/make_PON_and_run_mutect2/final_scripts_runs_DSL2/results_mutect2_filtering"


process getFilteredVariants_and_annotate {

	tag "FilterMutectCalls $bam_id"
        publishDir "$params.outdirB/filtered/", mode:'copy'

	input:
		path ${bam_id}-T_${bam_id}-N.filtered.vcf.gz	
		tuple val(bam_id) , file(bams)

	output:


	gatk IndexFeatureFile \
     		--input ${bam_id}-T_${bam_id}-N.filtered.vcf.gz

	gatk SelectVariants \
     		-R $refdir/Homo_sapiens_assembly38.fasta \
     		-V ${bam_id}-T_${bam_id}-N.filtered.vcf.gz \
     		--exclude-filtered true \
     		-O ${bam_id}-T_${bam_id}-N.filtered_only.vcf.gz


	}
