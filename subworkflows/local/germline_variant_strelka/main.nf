// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { STRELKA_GERMLINE                                              } from '../../../modules/nf-core/strelka/germline/main'
include { GATK4_MERGEVCFS as MERGE_STRELKA                              } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_INDEL                      } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_SNV                        } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_INDEL_GERMLINEVARIANT  } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_SNV_GERMLINEVARIANT    } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_INDEL                  } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_SNV                    } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                                                       } from '../../../modules/icgc-argo-workflows/cleanup/main'
workflow GERMLINE_VARIANT_STRELKA {

    take:
    // TODO nf-core: edit input (take) channels
        cram
        dict
        fasta
        fasta_fai
        analysis_json
        versions

    main:

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(versions)

    //Run Caller
    //Note Srelka produces two VCF types : variants and genome
    //We only want variants
    STRELKA_GERMLINE(cram,fasta, fasta_fai)
    ch_versions = ch_versions.mix(STRELKA_GERMLINE.out.versions)

    //Split outputs based on intervals VCF & TBI
    STRELKA_GERMLINE.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{strelka_vcf}

    STRELKA_GERMLINE.out.vcf_tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{strelka_tbi}

    //If intervals collect and merge into single
    MERGE_STRELKA(
        strelka_vcf.intervals
            .map{ meta, vcf ->
                new_meta = [
                                id : meta.id,
                                experimentalStrategy : meta.experimentalStrategy,
                                genomeBuild : meta.genomeBuild,
                                tumourNormalDesignation : meta.tumourNormalDesignation,
                                sampleType : meta.sampleType ,
                                gender : meta.gender,
                                study_id : meta.study_id,
                            ]

                [groupKey(new_meta, meta.num_intervals), vcf]
            }.groupTuple(),
        dict.map{ it -> [[id:it[0].baseName], it]}
    )
    ch_versions = ch_versions.mix(MERGE_STRELKA.out.versions)

    //Collect merged intervals VCF or single VCF
    merged_strelka_vcf = Channel.empty().mix(
                    MERGE_STRELKA.out.vcf,
                    strelka_vcf.no_intervals)
                .map{ meta, vcf ->
                    [[
                        id : meta.id,
                        experimentalStrategy : meta.experimentalStrategy,
                        genomeBuild : meta.genomeBuild,
                        tumourNormalDesignation : meta.tumourNormalDesignation,
                        sampleType : meta.sampleType ,
                        gender : meta.gender,
                        study_id : meta.study_id,
                        tool:  "strelka"
                    ],vcf]
                }

    //Same as above except for TBI
    merged_strelka_tbi = Channel.empty().mix(
                        MERGE_STRELKA.out.tbi,
                        strelka_tbi.no_intervals)
                    .map{ meta, tbi ->
                        [[
                            id : meta.id,
                            experimentalStrategy : meta.experimentalStrategy,
                            genomeBuild : meta.genomeBuild,
                            tumourNormalDesignation : meta.tumourNormalDesignation,
                            sampleType : meta.sampleType ,
                            gender : meta.gender,
                            study_id : meta.study_id,
                            tool:  "strelka"
                        ], tbi]
                    }

    // Recombine VCF and TBI to make a single payload
    select_vcf=merged_strelka_vcf
    .combine(merged_strelka_tbi)
    .map{ meta, vcf, metaB, tbi ->
        [[
            id: meta.id,
            experimentalStrategy : meta.experimentalStrategy,
            genomeBuild : meta.genomeBuild,
            tumourNormalDesignation : meta.tumourNormalDesignation,
            sampleType : meta.sampleType ,
            gender : meta.gender,
            study_id : meta.study_id,
            num_intervals:  meta.num_intervals,
            tool : meta.tool
        ],vcf,tbi,[]]
    }

    //Filter for Indels
    SELECT_VCF_INDEL(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_INDEL.out.versions)


    //Filter for SNVs
    SELECT_VCF_SNV(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_SNV.out.versions)
                    
    //Manipulate for payload ingestion
    ch_payload_vcf_indel=SELECT_VCF_INDEL.out.vcf.combine(SELECT_VCF_INDEL.out.tbi).combine(analysis_json)
            .map { metaA,vcf,metaB,tbi,analysis_json->
            [
                [
                    id : metaA.id,
                    experimentalStrategy : metaA.experimentalStrategy,
                    genomeBuild : metaA.genomeBuild,
                    tumourNormalDesignation : metaA.tumourNormalDesignation,
                    sampleType : metaA.sampleType ,
                    gender : metaA.gender,
                    study_id : metaA.study_id,
                    dataType : "InDel",
                    tool : metaA.tool
                ]
                ,[vcf, tbi],analysis_json]
            }

    ch_payload_vcf_snv=SELECT_VCF_SNV.out.vcf.combine(SELECT_VCF_SNV.out.tbi).combine(analysis_json)
            .map { metaA,vcf,metaB,tbi,analysis_json->
            [
                [
                    id : metaA.id,
                    experimentalStrategy : metaA.experimentalStrategy,
                    genomeBuild : metaA.genomeBuild,
                    tumourNormalDesignation : metaA.tumourNormalDesignation,
                    sampleType : metaA.sampleType ,
                    gender : metaA.gender,
                    study_id : metaA.study_id,
                    dataType : "SNV",
                    tool : metaA.tool
                ]
                ,[vcf, tbi],analysis_json]
            }

    //Generate payload
    PAYLOAD_VCF_INDEL_GERMLINEVARIANT(
        ch_payload_vcf_indel,
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        false
    )

    PAYLOAD_VCF_SNV_GERMLINEVARIANT(
        ch_payload_vcf_snv,
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        false
    )

    //Gather temporary files
    ch_cleanup=MERGE_STRELKA.out.vcf.map{meta,vcf -> vcf}
        .mix(STRELKA_GERMLINE.out.vcf.map{meta,vcf -> [vcf]})
        .mix(SELECT_VCF_INDEL.out.vcf.map{meta,vcf -> [vcf]})
        .mix(SELECT_VCF_SNV.out.vcf.map{meta,vcf -> [vcf]})
        .mix(PAYLOAD_VCF_INDEL_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
        .mix(PAYLOAD_VCF_SNV_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
        .collect()


    //If Local is true, will be published into "output_dir" directory
    if (params.local==false){
        //Upload variants
        SONG_SCORE_UP_VCF_INDEL(PAYLOAD_VCF_INDEL_GERMLINEVARIANT)
        SONG_SCORE_UP_VCF_SNV(PAYLOAD_VCF_SNV_GERMLINEVARIANT)
        if (params.cleanup){
            CLEANUP(
                ch_cleanup.collect(),
                SONG_SCORE_UP_VCF_INDEL.out.analysis_id
                .combine(SONG_SCORE_UP_VCF_SNV.out.analysis_id).collect()
            )
        }
    } else {
     if (params.cleanup){
        CLEANUP(ch_cleanup.collect(),
        PAYLOAD_VCF_INDEL_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]}
        .combine(PAYLOAD_VCF_SNV_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
        )
        }       
    }


    emit:
    versions = ch_versions
}

