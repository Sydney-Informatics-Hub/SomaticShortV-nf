#!/usr/bin/env nextflow
nextflow.enable.dsl = 2




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



include { mutect2 } from './mutect2'
include { GatherVcfs_step } from './GatherVcfs_step'

include { MergeMutectStats } from './MergeMutectStats'
include { LearnReadOrientationModel } from './LearnReadOrientationModel'

include { GetPileupSummaries_T; GetPileupSummaries_N } from './GetPileupSummaries'

include{ CalculateContamination } from './CalculateContamination'

include { FilterMutectCalls } from './FilterMutectCalls'

include { getFilteredVariants_and_annotate } from './getFilteredVariants_and_annotate'


workflow {

	mutect2(params.ponvcf,params.ponvcf_index,bam_pair_ch,intervalList,base_path,refdir)

	GatherVcfs_step(mutect2.out[0].collect(),bam_pair_ch,base_path,refdir)

	MergeMutectStats(bam_pair_ch,GatherVcfs_step.out[0].collect(),base_path,refdir)
	
	LearnReadOrientationModel(MergeMutectStats.out[1].collect(),bam_pair_ch,base_path,refdir)

	GetPileupSummaries_T(params.common_biallelic_path,params.common_biallelic_idx_path,bam_pair_ch,LearnReadOrientationModel.out.collect(),base_path,refdir)
	GetPileupSummaries_N(params.common_biallelic_path,params.common_biallelic_idx_path,bam_pair_ch,LearnReadOrientationModel.out.collect(),base_path,refdir)


	CalculateContamination(bam_pair_ch,GetPileupSummaries_T.out.collect(),GetPileupSummaries_N.out.collect(),base_path,refdir)

	FilterMutectCalls(bam_pair_ch,CalculateContamination.out[0].collect(),params.outdir,base_path,refdir)	

	getFilteredVariants_and_annotate(bam_pair_ch,FilterMutectCalls.out.collect(),params.outdir,base_path,refdir)


	}


