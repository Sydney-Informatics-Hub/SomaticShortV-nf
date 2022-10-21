#!/usr/bin/env nextflow
  


command_path="/scratch/wz54/npd561/installations/snpEff/snpEff"


process getFilteredVariants_and_annotate {

        tag "getFilteredVariants_and_annotate $bam_id"
        publishDir "$params.outdir/Mutect2_filtered/", mode:'copy'


        input:
                tuple val(bam_id) , file(bams)
                path pair_filtered_vcfs_allvariants

                path outdir

                path base_path
                path refdir


        output:
                path ("${bam_id}-T_${bam_id}-N.filtered_only.vcf.gz")
                path("${bam_id}-T_${bam_id}-N.filtered_only.ann.vcf")

        shell:

        '''
        gatk IndexFeatureFile \
                --input !{params.outdir}/Mutect2_filtered/!{bam_id}-T_!{bam_id}-N.filtered.vcf.gz

        gatk SelectVariants \
                -R !{refdir}/Homo_sapiens_assembly38.fasta \
                -V !{params.outdir}/Mutect2_filtered/!{bam_id}-T_!{bam_id}-N.filtered.vcf.gz \
                --exclude-filtered true \
                -O !{bam_id}-T_!{bam_id}-N.filtered_only.vcf.gz


        java -Xmx8g -jar !{command_path}/snpEff.jar -v \
                -o gatk \
                -stats !{bam_id}-T_!{bam_id}-N.filtered_only.ann.html \
                GRCh38.86 \
                !{bam_id}-T_!{bam_id}-N.filtered_only.vcf.gz > !{bam_id}-T_!{bam_id}-N.filtered_only.ann.vcf

        '''

        }
