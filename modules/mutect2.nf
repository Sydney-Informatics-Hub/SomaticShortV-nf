#!/usr/bin/env nextflow


refdir='/scratch/wz54/gs5517/sarek_testing/Reference/v0'
base_path='/scratch/wz54/npd561/PIPE-2629_thyroid_carcinoma/nextflow_pipelines/nfcore/sarek'

process mutect2 {


        tag "$bam_id $splitIntervalNumber"

        publishDir "$params.outdirA/", mode:'copy'

        input:
                path pon_vcf
                path pon_vcf_index
                tuple val(bam_id) , file(bams)

                each splitIntervalNumber

        output:
                path ("${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz")
                path ("${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz.stats")
                path ("${bam_id}-T_${bam_id}-N.f1r2.${splitIntervalNumber}.tar.gz")


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

