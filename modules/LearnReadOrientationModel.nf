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
  tag "LearnReadOrientationModel $bam_id"
  publishDir "$params.outdir/", mode:'copy'

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
	  path ("${bam_id}-T_${bam_id}-N.read-orientation-model.tar.gz")
	
	
	
	// this is the code block 
	shell:
    '''
      # Change this when running the complete ' mutect2' pipeline - followed by these filtering steps 
      ls !{base_path}/results_mutect2/!{bam_id}*f1r2.*.tar.gz > !{bam_id}-T_!{bam_id}-N.unfiltered_f1r2.args
  
      gatk LearnReadOrientationModel --java-options "-Xmx58g" \
        --input !{bam_id}-T_!{bam_id}-N.unfiltered_f1r2.args \
        -O !{bam_id}-T_!{bam_id}-N.read-orientation-model.tar.gz
    '''

}
