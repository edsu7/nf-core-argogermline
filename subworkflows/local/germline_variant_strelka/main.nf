// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { STRELKA_GERMLINE                              } from '../../../modules/nf-core/strelka/germline/main'
include { GATK4_MERGEVCFS as MERGE_STRELKA              } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { PAYLOAD_GERMLINEVARIANT                       } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD                             } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                                       } from '../../../modules/icgc-argo-workflows/cleanup/main'
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
                                id:             meta.id,
                                num_intervals:  meta.num_intervals
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
                        id:             meta.id,
                        num_intervals:  meta.num_intervals,
                        variantcaller:  "strelka"
                    ],vcf]
                }

    //Same as above except for TBI
    merged_strelka_tbi = Channel.empty().mix(
                        MERGE_STRELKA.out.tbi,
                        strelka_tbi.no_intervals)
                    .map{ meta, tbi ->
                        [[
                            id:             meta.id,
                            num_intervals:  meta.num_intervals,
                            variantcaller:  "deepvariant"
                        ], tbi]
                    }
                    
    //Manipulate for payload ingestion
    ch_payload=merged_strelka_vcf.combine(merged_strelka_tbi).combine(analysis_json)
            .map { metaA,vcf,metaB,tbi,analysis_json->
            [
                [ id : metaA.id,
                  study_id : params.study_id,
                  tool : "strelka"
                ]
                ,[vcf, tbi],analysis_json]
            }

    //Generate payload
    PAYLOAD_GERMLINEVARIANT(
        ch_payload,
        "",
        "",
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        "strelka"
    )

    //Gather temporary files
        ch_cleanup=MERGE_STRELKA.out.vcf.map{meta,vcf -> vcf}
        .mix(STRELKA_GERMLINE.out.vcf.map{meta,vcf -> [vcf]})
        .mix(PAYLOAD_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
        .collect()

    //If Local is true, will be published into "output_dir" directory
    if (params.local==false){
        //Upload variants
        //SONG_SCORE_UPLOAD(PAYLOAD_GERMLINEVARIANT.out.payload_files)
        if (params.cleanup){
            CLEANUP(ch_cleanup.collect(),SONG_SCORE_UPLOAD.analysis_id)
        }
    } else {
     if (params.cleanup){
        CLEANUP(ch_cleanup.collect(),PAYLOAD_GERMLINEVARIANT.out.payload_files)
        }       
    }

    emit:
    // TODO nf-core: edit emitted channels
    versions = ch_versions
}

