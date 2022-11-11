#!/usr/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2

// Define the process
/// This process runs mutect2 on a tumor/normal sample-pair 

process mutect2 {

        // Unhash container command below and edit name of container
	// if using Docker/Singularity containers
        //container "${params.container}
        
        // where to publish the outputs
        tag "$bam_id $splitIntervalNumber"
        publishDir "$params.outdirA/", mode:'copy'
        
        
        // See: https://www.nextflow.io/docs/latest/process.html#inputs
	/// each input needs to be placed on a new line
        input:
                path pon_vcf
                path pon_vcf_index
                tuple val(bam_id) , file(bams)

                each splitIntervalNumber

        // See: https://www.nextflow.io/docs/latest/process.html#outputs
	// each new output needs to be placed on a new line
        output:
                path ("${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz")
                path ("${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz.stats")
                path ("${bam_id}-T_${bam_id}-N.f1r2.${splitIntervalNumber}.tar.gz")

        // this is the code block 
        script:

        """

        gatk Mutect2 \
             -R $refdir/Homo_sapiens_assembly38.fasta \
             -I ${bams[1]} \
             -I ${bams[0]} \
             -normal ${bam_id}-N \
             --panel-of-normals ${pon_vcf} \
             --germline-resource $refdir/gnomAD.r2.1.1.GRCh38.PASS.AC.AF.only.vcf.gz \
             --f1r2-tar-gz ${bam_id}-T_${bam_id}-N.f1r2.${splitIntervalNumber}.tar.gz \
             -XL chrM \
             -L "$base_path/Somatic-ShortV/nextflow/make_PON_and_run_mutect2/final_scripts_runs_DSL2/100M_primary_interval_${splitIntervalNumber}.list" \
             -O ${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz
        """




}

