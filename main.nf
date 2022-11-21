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
          SOMATIC SHORT V - NF 
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
  Usage:  nextflow run main.nf --ref reference.fasta

  Required Arguments:
	//--input		  Full path and name of sample input file (tsv format).
	--ref			  Full path and name of reference genome (fasta format).
	
  Optional Arguments:
  --outDir    Specify name of results directory. 
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
	
	
	// Set PATH pointing to the input 'bam' file-pairs 
	params.bams = "$base_dir/Preprocessing/*/Recalibrated/*-{N,T}.recal.bam"
	// bam pair channel
	bam_pair_ch=Channel.fromFilePairs( params.bams )
	
	
	



	# Run the processes 
	mutect2(params.ponvcf,params.ponvcf+'.tbi',bam_pair_ch,intervalList,base_path,refdir_path)

	GatherVcfs_step(mutect2.out[0].collect(),bam_pair_ch)

	MergeMutectStats(bam_pair_ch,GatherVcfs_step.out[0].collect(),base_path)

	LearnReadOrientationModel(MergeMutectStats.out[1].collect(),bam_pair_ch,base_path)

	GetPileupSummaries_T(params.common_biallelic_path,params.common_biallelic_path+'.tbi', bam_pair_ch,LearnReadOrientationModel.out.collect())

	GetPileupSummaries_N(params.common_biallelic_path,params.common_biallelic_path+'.tbi',bam_pair_ch,LearnReadOrientationModel.out.collect())

	CalculateContamination(bam_pair_ch,GetPileupSummaries_T.out.collect(),GetPileupSummaries_N.out.collect())

	FilterMutectCalls(bam_pair_ch,CalculateContamination.out[0].collect(),params.outdir)	

	getFilteredVariants_and_annotate(bam_pair_ch,FilterMutectCalls.out.collect(),params.outdir,refdir_path)


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
