#!/bin/env nextflow 

// Enable DSL-2 syntax
nextflow.enable.dsl=2

// Define the process
/// This process Gathers multiple VCF files from a scatter operation into a single VCF file. 
/// Input files must be supplied in genomic order and must not have events at overlapping positions.

process GatherVcfs_step {
  cpus "${params.cpus}"
	debug = true

	// Unhash container command below and edit name of container
	// if using Docker/Singularity containers
        //container "${params.container}

	// where to publish the outputs
        tag "GatherVcfs_step $bam_id"
        publishDir "$params.outdir/Mutect2_unfiltered", mode:'copy'

	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	/// each input needs to be placed on a new line
  input:
    path ('*')
    tuple val(bam_id) , file(bams)
    path base_path
	  path refdir

	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	// each new output needs to be placed on a new line
	output:
	  path ("${bam_id}-T_${bam_id}-N.unfiltered.vcf.gz")
    path ("${bam_id}-T_${bam_id}-N.unfiltered.vcf.gz.tbi")
	
	
	
	// this is the code block 
	script:
	"""
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
