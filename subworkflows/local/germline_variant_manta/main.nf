// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { MANTA_GERMLINE                                                } from '../../../modules/nf-core/manta/germline/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_DIPLOID                        } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_INDEL                      } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_SNV                        } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_INDEL_GERMLINEVARIANT  } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_SNV_GERMLINEVARIANT    } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_INDEL                  } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_SNV                    } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                                                       } from '../../../modules/icgc-argo-workflows/cleanup/main'

workflow GERMLINE_VARIANT_MANTA {

    take:
        cram                     // channel: [mandatory] [meta, cram, crai, interval.bed.gz, interval.bed.gz.tbi]
        dict                     // channel: [optional]
        fasta                    // channel: [mandatory]
        fasta_fai                // channel: [mandatory]
        analysis_json
        versions

    main:

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(versions)

    //Run caller
    MANTA_GERMLINE(
        cram,
        fasta.map{ it -> [[id:it[0].baseName], it]}, 
        fasta_fai.map{ it -> [[id:it[0].baseName], it]}
        )
    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions.first())

    //Split according to VCF and TBI by intervals
    MANTA_GERMLINE.out.diploid_sv_vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_diploid_sv_vcf}
    MANTA_GERMLINE.out.diploid_sv_vcf_tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{manta_diploid_sv_vcf_tbi}

    // Only when using intervals
    MERGE_MANTA_DIPLOID(
        manta_diploid_sv_vcf.intervals
            .map{ meta, vcf ->

                [groupKey([
                            id: meta.id,
                            experimentalStrategy : meta.experimentalStrategy,
                            genomeBuild : meta.genomeBuild,
                            tumourNormalDesignation : meta.tumourNormalDesignation,
                            sampleType : meta.sampleType ,
                            gender : meta.gender,
                            study_id : meta.study_id,
                            num_intervals:  meta.num_intervals
                        ],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict.map{ it -> [[id:it[0].baseName], it]})
    ch_versions = ch_versions.mix(MERGE_MANTA_DIPLOID.out.versions)

    // Mix output channels for "no intervals" and "with intervals" results
    // Only diploid SV should get annotated
    collected_manta_diploid_sv_vcf = Channel.empty().mix(
                    MERGE_MANTA_DIPLOID.out.vcf,
                    manta_diploid_sv_vcf.no_intervals)
                .map{ meta, vcf ->
                    [[
                        id: meta.id,
                        experimentalStrategy : meta.experimentalStrategy,
                        genomeBuild : meta.genomeBuild,
                        tumourNormalDesignation : meta.tumourNormalDesignation,
                        sampleType : meta.sampleType ,
                        gender : meta.gender,
                        study_id : meta.study_id,
                        num_intervals:  meta.num_intervals,
                        tool:  "manta"],
                    vcf]
                }
    collected_manta_diploid_sv_vcf_tbi = Channel.empty().mix(
                    MERGE_MANTA_DIPLOID.out.tbi,
                    manta_diploid_sv_vcf_tbi.no_intervals)
                .map{ meta, tbi ->
                    [[
                        id: meta.id,
                        experimentalStrategy : meta.experimentalStrategy,
                        genomeBuild : meta.genomeBuild,
                        tumourNormalDesignation : meta.tumourNormalDesignation,
                        sampleType : meta.sampleType ,
                        gender : meta.gender,
                        study_id : meta.study_id,
                        num_intervals:  meta.num_intervals,
                        tool:  "manta"],
                    tbi]
                }

    // Recombine VCF and TBI to make a single payload
    select_vcf=collected_manta_diploid_sv_vcf
    .combine(collected_manta_diploid_sv_vcf_tbi)
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
    select_vcf.view()
    //Filter for Indels
    SELECT_VCF_INDEL(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_INDEL.out.versions)


    //Filter for SNVs
    SELECT_VCF_SNV(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_SNV.out.versions)

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
    ch_payload_vcf_indel.view()
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
    ch_cleanup=MANTA_GERMLINE.out.diploid_sv_vcf.map{meta,vcf -> [vcf]}
    .mix(MERGE_MANTA_DIPLOID.out.vcf.map{meta,vcf -> [vcf]})
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

