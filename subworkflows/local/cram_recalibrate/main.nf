// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { GATK4_BASERECALIBRATOR                                } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_GATHERBQSRREPORTS                               } from '../../../modules/nf-core/gatk4/gatherbqsrreports/main'
include { GATK4_APPLYBQSR                                       } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX                                        } from '../../../modules/nf-core/samtools/index/main'
include { CLEANUP                                               } from '../../../modules/icgc-argo-workflows/cleanup/main'
workflow CRAM_RECALIBRATE {

    take:
        cram
        dict
        fasta
        fasta_fai
        intervals
        known_sites
        known_sites_tbi

    main:

        ch_versions = Channel.empty()

        //Combine CRAM and Interval channels
        cram_intervals = cram.combine(intervals)
        .map{ meta, cram, crai, intervals, num_intervals ->
                intervals_new = num_intervals == 0 ? [] : intervals
                [[
                    experimentalStrategy : meta.experimentalStrategy,
                    genomeBuild : meta.genomeBuild ,
                    tumourNormalDesignation : meta.tumourNormalDesignation,
                    sampleType : meta.sampleType,
                    gender : meta.gender,
                    id : meta.id,
                    num_intervals:  num_intervals
                ],
                cram, crai, intervals_new]
            }

        //Recalibrate
        GATK4_BASERECALIBRATOR(cram_intervals, fasta, fasta_fai, dict, known_sites, known_sites_tbi)
        ch_versions=ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

        ch_cram_applybqsr = cram_intervals.combine(GATK4_BASERECALIBRATOR.out.table)               

        //Combine recalibrate table output with cram and intervals
        ch_cram_applybqsr = cram_intervals.combine(GATK4_BASERECALIBRATOR.out.table)
            .map{ meta, cram, crai, intervals, num_intervals,recal ->
                 [[
                    experimentalStrategy : meta.experimentalStrategy,
                    genomeBuild : meta.genomeBuild ,
                    tumourNormalDesignation : meta.tumourNormalDesignation,
                    sampleType : meta.sampleType,
                    gender : meta.gender,
                    id : meta.id,
                 ],
                 cram,crai,recal,intervals]
             }

        //Apply base recalibration based on recal.table
        GATK4_APPLYBQSR(ch_cram_applybqsr,fasta,fasta_fai,dict)
        ch_versions=ch_versions.mix(GATK4_APPLYBQSR.out.versions)

        //Index recalibrated CRAM
        SAMTOOLS_INDEX(GATK4_APPLYBQSR.out.cram)
        ch_versions=ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        //Combine previous yml into singular, this avoids conflicts as previous ymls are cleaned up
        ch_versions=ch_versions.unique().collectFile(name: 'recalibrate_versions.yml')

        //Combine cram and crai into single channel
        cram_crai = GATK4_APPLYBQSR.out.cram.combine(SAMTOOLS_INDEX.out.crai)
            .map{ metaA, cram, metaB, crai ->
            [[
                experimentalStrategy : metaA.experimentalStrategy,
                genomeBuild : metaA.genomeBuild ,
                tumourNormalDesignation : metaA.tumourNormalDesignation,
                sampleType : metaA.sampleType,
                gender : metaA.gender,
                id : metaA.id,
            ],
            cram,crai]
        }

        //Collect temporary files and remove 
        ch_cleanup=cram.map{meta,cram,crai -> [cram]}.collect()
        .mix(GATK4_BASERECALIBRATOR.out.table.map{meta,table -> [table]}.collect()).collect()


        if (params.cleanup){
            CLEANUP(ch_cleanup.collect(),cram_crai)
        }   


    emit:
        cram = cram_crai
        versions = ch_versions 
}

