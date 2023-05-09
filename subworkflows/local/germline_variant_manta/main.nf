// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { MANTA_GERMLINE                         } from '../../../modules/nf-core/manta/germline/main'
include { GATK4_MERGEVCFS as MERGE_MANTA_DIPLOID } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { PAYLOAD_GERMLINEVARIANT                } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD                      } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                                } from '../../../modules/icgc-argo-workflows/cleanup/main'

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
    MANTA_GERMLINE(cram, fasta, fasta_fai)
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
                            id:             meta.id,
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
                        id:             meta.id,
                        num_intervals:  meta.num_intervals,
                        variantcaller:  "manta"],
                    vcf]
                }
    collected_manta_diploid_sv_vcf_tbi = Channel.empty().mix(
                    MERGE_MANTA_DIPLOID.out.tbi,
                    manta_diploid_sv_vcf_tbi.no_intervals)
                .map{ meta, tbi ->
                    [[
                        id:             meta.id,
                        num_intervals:  meta.num_intervals,
                        variantcaller:  "manta"],
                    tbi]
                }

    //Manipulate for payload ingestion
    ch_payload=collected_manta_diploid_sv_vcf
    .combine(collected_manta_diploid_sv_vcf_tbi)
    .combine(analysis_json)
    .map {metaA,vcf,metaB,tbi,analysis_json ->
    [
        [ id : metaA.id,
            study_id : params.study_id,
            tool : "manta"
        ],
        [vcf, tbi], analysis_json]
    }

    //Generate payload
    PAYLOAD_GERMLINEVARIANT(
        ch_payload,
        "",
        "",
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        "manta"
    )

    //Gather temporary files
    ch_cleanup=MANTA_GERMLINE.out.diploid_sv_vcf.map{meta,vcf -> [vcf]}
    .mix(MERGE_MANTA_DIPLOID.out.vcf.map{meta,vcf -> [vcf]})
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
    versions = ch_versions                     // channel: [ versions.yml ]
}

