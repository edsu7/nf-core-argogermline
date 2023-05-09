// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { DEEPVARIANT                               } from '../../../modules/nf-core/deepvariant/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_GVCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF  } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_GVCF  } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_VCF   } from '../../../modules/nf-core/tabix/tabix/main'
include { PAYLOAD_GERMLINEVARIANT                   } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD                         } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                                   } from '../../../modules/icgc-argo-workflows/cleanup/main'

workflow GERMLINE_VARIANT_DEEPVARIANT {

    take:
        cram                     // channel: [mandatory] [meta, cram, crai, interval]
        dict                     // channel: [optional]
        fasta                    // channel: [mandatory]
        fasta_fai                // channel: [mandatory]
        analysis_json
        versions

    main:

    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(versions)

    //Run caller
    DEEPVARIANT(
        cram,     // [meta, cram , crai]
        fasta,    // [meta , fasta]
        fasta_fai // [meta , fasta]
        )
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())


    //Separate based on intervals
    DEEPVARIANT.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }.set{deepvariant_vcf_out}

    // Only when no intervals, compress and index
    TABIX_VC_DEEPVARIANT_VCF(deepvariant_vcf_out.no_intervals)
    ch_versions = ch_versions.mix(TABIX_VC_DEEPVARIANT_VCF.out.versions)

    // When intervals, merge, compres and index
    MERGE_DEEPVARIANT_VCF(
        deepvariant_vcf_out.intervals
            .map{ meta, vcf ->

                new_meta = [
                            id:             meta.id,
                            num_intervals:  meta.num_intervals
                        ]
                [groupKey(new_meta, meta.num_intervals), vcf]
            }.groupTuple(),
        dict.map{ it -> [[id:it[0].baseName], it]})
    ch_versions = ch_versions.mix(MERGE_DEEPVARIANT_VCF.out.versions)

    //Combine single and multi interval VCF channels
    deepvariant_vcf = Channel.empty().mix(
                        MERGE_DEEPVARIANT_VCF.out.vcf,
                        deepvariant_vcf_out.no_intervals)
                    .map{ meta, vcf ->
                        [[
                            id:             meta.id,
                            num_intervals:  meta.num_intervals,
                            variantcaller:  "deepvariant"
                        ], vcf]
                    }
    //Combine single and multi interval TBI channels
    deepvariant_tbi = Channel.empty().mix(
                        MERGE_DEEPVARIANT_VCF.out.tbi,
                        TABIX_VC_DEEPVARIANT_VCF.out.tbi)
                    .map{ meta, tbi ->
                        [[
                            id:             meta.id,
                            num_intervals:  meta.num_intervals,
                            variantcaller:  "deepvariant"
                        ], tbi]
                    }

    //Manipulate for payload ingestion
    ch_payload=deepvariant_vcf.combine(deepvariant_tbi).combine(analysis_json)
            .map { metaA,vcf,metaB,tbi,analysis_json->
            [
                [ id : metaA.id,
                  study_id : params.study_id,
                  tool : "deepvariant"
                ]
                ,[vcf, tbi],analysis_json]
            }

    //Generate payload
    PAYLOAD_GERMLINEVARIANT(
        ch_payload,
        "",
        "",
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        "deepvariant"
    )

    //Gather temporary files
    ch_cleanup=DEEPVARIANT.out.vcf.map{meta,vcf -> [vcf]}
    .mix(TABIX_VC_DEEPVARIANT_VCF.out.tbi.map{meta,tbi -> [tbi]})
    .mix(MERGE_DEEPVARIANT_VCF.out.vcf.map{meta,vcf -> [vcf]})
    .mix(PAYLOAD_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]}).collect()

    //If Local is true, will be published into "output_dir" directory
    if (params.local==false){
        //Upload variants
        //SONG_SCORE_UPLOAD(PAYLOAD_GERMLINEVARIANT.out.payload_files)
        if (params.cleanup){
            CLEANUP(ch_cleanup.collect(),SONG_SCORE_UPLOAD.out.analysis_id)
        }
    } else {
     if (params.cleanup){
        CLEANUP(ch_cleanup.collect(),PAYLOAD_GERMLINEVARIANT.out.payload_files)
        }       
    }


    emit:
    versions = ch_versions
}

