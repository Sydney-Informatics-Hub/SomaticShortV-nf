#!/usr/bin/env nextflow


process FilterMutectCalls {

        tag "FilterMutectCalls $bam_id"
        publishDir "$params.outdir/filtered/", mode:'copy'


        input :
		tuple val(bam_id) , file(bams)
                path pair_contaminationTable

                path outdir

        output :
                path ("${bam_id}-T_${bam_id}-N.filtered.vcf.gz")

        shell:
        '''
        gatk --java-options "-Xmx16g -Xms16g" \
                FilterMutectCalls \
                --reference !{refdir}/Homo_sapiens_assembly38.fasta \
                -V !{params.outdirA}/Mutect2/!{bam_id}-T_!{bam_id}-N.unfiltered.vcf.gz \
                --stats !{params.outdir}/!{bam_id}-T_!{bam_id}-N.unfiltered.stats \
                --tumor-segmentation !{params.outdir}/!{bam_id}-T_segments.table \
                --contamination-table  !{params.outdir}/!{bam_id}-T_!{bam_id}-N_contamination.table \
                --ob-priors !{params.outdir}/!{bam_id}-T_!{bam_id}-N.read-orientation-model.tar.gz \
                -O !{bam_id}-T_!{bam_id}-N.filtered.vcf.gz



        '''


        }

