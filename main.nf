#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2



// =================================================================
// main.nf is the pipeline script for a nextflow pipeline
// Should contain the following sections:
	// Import subworkflows
	// Log info function
	// Help function 
	// Main workflow structure
	// Some tests to check input data, essential arguments

// Examples are included for each section. Remove them and replace
// with project-specific code. For more information on nextflow see:
// https://www.nextflow.io/docs/latest/index.html and the SIH Nextflow
// upskilling repo @ INSERT REPO PATH 
//
// ===================================================================

// Import subworkflows to be run in the workflow
// Each of these is a separate .nf script saved in ./modules/

  
include { mutect2 } from './modules/mutect2'
include { GatherVcfs_step } from './modules/GatherVcfs_step'
include { MergeMutectStats } from './modules/MergeMutectStats'
include { LearnReadOrientationModel } from './modules/LearnReadOrientationModel'
include { GetPileupSummaries_T; GetPileupSummaries_N } from './modules/GetPileupSummaries'
include{ CalculateContamination } from './modules/CalculateContamination'
include { FilterMutectCalls } from './modules/FilterMutectCalls'
include { getFilteredVariants_and_annotate } from './modules/getFilteredVariants_and_annotate'
  


/// Print a header for your pipeline 

log.info """\

      ============================
      ============================
          GERMLINE SHORT V - NF 
      ============================
      ============================

 -._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _  
    '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '.` :    
  '.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.:    
  : '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.  
  '   '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.  
         `-..,..-'       `-..,..-'       `-..,..-'       `       


             ~~~~ Version: 1.0 ~~~~
 

 Created by the Sydney Informatics Hub, University of Sydney

 Find documentation and more info @ GITHUB REPO DOT COM

 Cite this pipeline @ INSERT DOI

 Log issues @ https://github.com/Sydney-Informatics-Hub/SomaticShortV-nf/issues

 All of the default parameters are set in `nextflow.config`

 =======================================================================================
Workflow run parameters 
=======================================================================================

input       : ${params.input}
outDir      : ${params.outDir}
workDir     : ${workflow.workDir}

=======================================================================================

 """

/// Help function 
// This is an example of how to set out the help function that 
// will be run if run command is incorrect (if set in workflow) 
// or missing/  

def helpMessage() {
    log.info"""
  Usage:  nextflow run <PATH TO REPO>/myPipeline-nf <args> --outDir

  Required Arguments:
	--outDir	Specify path to output directory

	--input		Specify full path and name of sample
			input file (tab separated).
    """.stripIndent()
}

/// Main workflow structure. 

workflow {

// Show help message if --help is run or if any required params are not 
// provided at runtime

        if ( params.help || params.input == false )
	{   
        // Invoke the help function above and exit
              helpMessage()
              exit 1

        // consider adding some extra contigencies here.
        // could validate path of all input files in list?
        // could validate indexes for input files exist?
        // could validate indexes for reference exist?
        // confirm with each tool, any requirements for their run?

        // if none of the above are a problem, then run the workflow
	} 
	else 
	{
	
  	// Define params and input channels 
	
	// 
	base_path=""
	
	// Set PATH pointing to the input 'bam' file-pairs 
	params.bams = "$base_path/Preprocessing/*/Recalibrated/*-{N,T}.recal.bam"
	// bam pair channel
	bam_pair_ch=Channel.fromFilePairs( params.bams )
	
	refdir="$base_path/Reference/v0"

	
	// TBD Add the files below to the above refdir
	
	params.common_biallelic_path="$base_path/allelic_references/small_exac_common_3.hg38.vcf.gz"
	params.common_biallelic_idx_path="$base_path/allelic_references/small_exac_common_3.hg38.vcf.gz.tbi"
	
	params.outdir="$base_path/results_mutect2"

	// PATH to PoN (created previously)
	params.ponvcf="$base_path/pon.vcf.gz"
	params.ponvcf_index="$base_path/pon.vcf.gz.tbi"

	// Intervals for 'Scatter-Gather'
	intervalList=['a','b','c','d','e','f','g','h','i','j','k','l','m','n']


	# Run the processes 
	mutect2(params.ponvcf,params.ponvcf_index,bam_pair_ch,intervalList,base_path,refdir)

	GatherVcfs_step(mutect2.out[0].collect(),bam_pair_ch,base_path,refdir)

	MergeMutectStats(bam_pair_ch,GatherVcfs_step.out[0].collect(),base_path,refdir)

	LearnReadOrientationModel(MergeMutectStats.out[1].collect(),bam_pair_ch,base_path,refdir)

	GetPileupSummaries_T(params.common_biallelic_path,params.common_biallelic_idx_path,bam_pair_ch,LearnReadOrientationModel.out.collect(),base_path,refdir)

	GetPileupSummaries_N(params.common_biallelic_path,params.common_biallelic_idx_path,bam_pair_ch,LearnReadOrientationModel.out.collect(),base_path,refdir)

	CalculateContamination(bam_pair_ch,GetPileupSummaries_T.out.collect(),GetPileupSummaries_N.out.collect(),base_path,refdir)

	FilterMutectCalls(bam_pair_ch,CalculateContamination.out[0].collect(),params.outdir,base_path,refdir)	

	getFilteredVariants_and_annotate(bam_pair_ch,FilterMutectCalls.out.collect(),params.outdir,base_path,refdir)


	}}

workflow.onComplete {
  summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
outDir      : ${params.outDir}

=======================================================================================
  """
  println summary

}
