// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { TIDDIT_SV                         } from '../../../modules/nf-core/tiddit/sv/main'
include { TABIX_BGZIPTABIX                  } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_INDEL                      } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_SNV                        } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_INDEL_GERMLINEVARIANT  } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_SNV_GERMLINEVARIANT    } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_INDEL                  } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_SNV                    } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                           } from '../../../modules/icgc-argo-workflows/cleanup/main'

workflow GERMLINE_VARIANT_TIDDIT {

    take:
        cram
        fasta
        bwa
        analysis_json
        versions

    main:
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(versions)

    //Run caller
    TIDDIT_SV(
        cram,  // [meta, cram , crai]
        fasta, // [meta , fasta]
        bwa    // [meta , bwa]
        )
    ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)

    //Compress and index
    TABIX_BGZIPTABIX(TIDDIT_SV.out.vcf)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    //Manipulate for payload ingestion
    TABIX_BGZIPTABIX.out.gz_tbi.combine(analysis_json)
            .map { meta,vcf,tbi,analysis_json->
            [
                [
                    id : meta.id,
                    experimentalStrategy : meta.experimentalStrategy,
                    genomeBuild : meta.genomeBuild,
                    tumourNormalDesignation : meta.tumourNormalDesignation,
                    sampleType : meta.sampleType ,
                    gender : meta.gender,
                    study_id : params.study_id,
                    tool : "tiddit"
                ]
                ,[vcf, tbi],analysis_json]
            }.set{ch_payload}

    // Recombine VCF and TBI to make a single payload
    select_vcf=TABIX_BGZIPTABIX.out.gz_tbi
    .map{ meta, vcf, tbi ->
        [[
            id: meta.id,
            experimentalStrategy : meta.experimentalStrategy,
            genomeBuild : meta.genomeBuild,
            tumourNormalDesignation : meta.tumourNormalDesignation,
            sampleType : meta.sampleType ,
            gender : meta.gender,
            study_id : meta.study_id,
            num_intervals:  meta.num_intervals,
            tool : "tiddit"
        ],vcf,tbi,[]]
    }

    //Filter for Indels
    SELECT_VCF_INDEL(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_INDEL.out.versions)


    //Filter for SNVs
    SELECT_VCF_SNV(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_SNV.out.versions)


    //Generate payload
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
    ch_cleanup=TIDDIT_SV.out.vcf.map{ meta,vcf -> vcf}
        .mix(TABIX_BGZIPTABIX.out.gz_tbi.map{metaB,vcf_gz,vcf_gz_tbi -> vcf_gz})
        .mix(SELECT_VCF_INDEL.out.vcf.map{meta,vcf -> [vcf]})
        .mix(SELECT_VCF_SNV.out.vcf.map{meta,vcf -> [vcf]})
        .mix(PAYLOAD_VCF_INDEL_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
        .mix(PAYLOAD_VCF_SNV_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})

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

