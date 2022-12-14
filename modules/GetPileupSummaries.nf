#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2

container "${params.gatk4__container}"

process GetPileupSummaries_T {


        tag "GetPileupSummaries $bam_id"
        publishDir "$params.outdirB/", mode:'copy'


        input:
                path common_biallelic
                path common_biallelic_idx

		tuple val(bam_id) , file(bams)

                path f1r2Args_T


        output:
                path ("${bam_id}-T.pileups.table")

        shell:
        '''
        gatk --java-options "-Xmx56g -XX:ParallelGCThreads=2" \
                GetPileupSummaries \
                -I !{bams[1]} \
                -V !{common_biallelic} \
                -L !{common_biallelic} \
                -O !{bam_id}-T.pileups.table


        '''




}


process GetPileupSummaries_N {


        tag "GetPileupSummaries $bam_id"
        publishDir "$params.outdirB/", mode:'copy'


        input:
                path common_biallelic
                path common_biallelic_idx

		tuple val(bam_id) , file(bams)

                path f1r2Args_N


        output:
                path ("${bam_id}-N.pileups.table")


        shell:
        '''

        gatk --java-options "-Xmx56g -XX:ParallelGCThreads=2" \
                GetPileupSummaries \
                -I !{bams[0]} \
                -V !{common_biallelic} \
                -L !{common_biallelic} \
                -O !{bam_id}-N.pileups.table


        '''

}

