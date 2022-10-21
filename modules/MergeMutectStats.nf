#!/bin/env nextflow 

// Enable DSL-2 syntax
nextflow.enable.dsl=2

// Define the process
/// This process combine the stats files across these scattered intervals 

process MergeMutectStats {
  cpus "${params.cpus}"
  debug = true

	// Unhash container command below and edit name of container
	// if using Docker/Singularity containers
        //container "${params.container}

	// where to publish the outputs
  tag "MergeMutectStats"
  publishDir "$params.outdir/", mode:'copy'

	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	// each input needs to be placed on a new line
  input:
    tuple val(bam_id) , file(bams)
    path ('*') 
    path base_path
		path refdir

	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	// each new output needs to be placed on a new line
	output:
		path ("${bam_id}-T_${bam_id}-N.unfiltered_stats.args")
    path ("${bam_id}-T_${bam_id}-N.unfiltered.stats")
	
	
  // this is the code block ; Notice the single quote for shell block (this enables both nextflow and bash variables to be inlcuded in the script)

	shell:
  '''

    ls !{params.outdir}/!{bam_id}*.stats   >!{bam_id}-T_!{bam_id}-N.unfiltered_stats.args
    gatk MergeMutectStats \
      --stats !{bam_id}-T_!{bam_id}-N.unfiltered_stats.args \
      -O !{bam_id}-T_!{bam_id}-N.unfiltered.stats
  '''
}
