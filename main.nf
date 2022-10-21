#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2


refdir='/scratch/wz54/gs5517/sarek_testing/Reference/v0'

// Could not set it to "" or "./" :-( - Some issue on the Nimbus instance - Will check this
base_path="/scratch/wz54/npd561/PIPE-2629_thyroid_carcinoma/nextflow_pipelines/nfcore/sarek/Somatic-ShortV/nextflow/make_PON_and_run_mutect2/final_scripts_production_withsnpEff/run_mutect2_and_filter_DSL2"

// TBD Add the files below to the above refdir
small_exac_common_path="/scratch/wz54/npd561/PIPE-2629_thyroid_carcinoma/nextflow_pipelines/nfcore/sarek/Somatic-ShortV/nextflow/make_PON_and_run_mutect2"
params.common_biallelic_path="$small_exac_common_path/small_exac_common_3.hg38.vcf.gz"
params.common_biallelic_idx_path="$small_exac_common_path/small_exac_common_3.hg38.vcf.gz.tbi"

// Set PATH pointing to the 'bam' files 
params.bams = "/scratch/er01/PIPE-2629-ThyroidCancer/nf_sarek/preprocess_*/Preprocessing/*/Recalibrated/*-{N,T}.recal.bam"
// bam pair channel
bam_pair_ch=Channel.fromFilePairs( params.bams )

params.outdir="$base_path/results_mutect2"

// PATH to PoN (created previously)
params.ponvcf="$base_path/pon.vcf.gz"
params.ponvcf_index="$base_path/pon.vcf.gz.tbi"

// Intervals for 'Scatter-Gather'
intervalList=['a','b','c','d','e','f','g','h','i','j','k','l','m','n']






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
// Each of these is a separate .nf script saved in modules/
// Add as many of these as you need. The example below will
// look for the process called process in modules/moduleName.nf
// Include { process } from './modules/moduleName'
  
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

      ============================================
      ============================================
      A DSL2 pipeline to identify somatic short-variants 
      ============================================
      ============================================

 -._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _  
    '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '.` :    
  '.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.:    
  : '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.  
  '   '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.  
         `-..,..-'       `-..,..-'       `-..,..-'       `       


             ~~~~ Version: ${params.version} ~~~~
 

 Created by the Sydney Informatics Hub, University of Sydney

 Find documentation and more info @ GITHUB REPO DOT COM

 Cite this pipeline @ INSERT DOI

 Log issues @ GITHUB REPO DOT COM

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

/// Main workflow structure. Include some input/runtime tests here.
// Make sure to comment what each step does for readability. 

workflow {

// Show help message if --help is run or if any required params are not 
// provided at runtime

        if ( params.help || params.input == false ){   
        // Invoke the help function above and exit
              helpMessage()
              exit 1

        // consider adding some extra contigencies here.
        // could validate path of all input files in list?
        // could validate indexes for input files exist?
        // could validate indexes for reference exist?
        // confirm with each tool, any requirements for their run?

// if none of the above are a problem, then run the workflow
	} else {
	
  // Define input channels 
  // cohort_ch = Channel.fromPath("${params.cohort}")
  // outDir_ch = Channel.fromPath("${params.outDir}")

// Run process 1 example
	processOne(cohort_ch, outDir_ch)
	
	// process 2 example 
	processTwo(processOne.out)

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
