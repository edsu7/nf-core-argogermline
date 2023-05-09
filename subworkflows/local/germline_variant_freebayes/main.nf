// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { FREEBAYES                                     } from '../../../modules/nf-core/freebayes/main'
include { GATK4_MERGEVCFS as MERGE_FREEBAYES            } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { TABIX_TABIX as TABIX_VC_FREEBAYES             } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_SORT                                 } from '../../../modules/nf-core/bcftools/sort/main'
include { PAYLOAD_GERMLINEVARIANT                       } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD                             } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                                       } from '../../../modules/icgc-argo-workflows/cleanup/main'

workflow GERMLINE_VARIANT_FREEBAYES {

    take:
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
    FREEBAYES(
        cram,
        fasta,
        fasta_fai,
        [], [], [])
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    BCFTOOLS_SORT(FREEBAYES.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    //Split outputs based on intervals VCF & TBI
    BCFTOOLS_SORT.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{bcftools_vcf_out}

    //Index if no intervals
    TABIX_VC_FREEBAYES(bcftools_vcf_out.no_intervals)


    MERGE_FREEBAYES(
        bcftools_vcf_out.intervals
            .map{ meta, vcf ->

                [groupKey([
                            id:             meta.id,
                            num_intervals:  meta.num_intervals
                        ],
                        meta.num_intervals),
                vcf]

            }.groupTuple(),
        dict.map{ it -> [[id:it[0].baseName], it]})
    ch_versions = ch_versions.mix(MERGE_FREEBAYES.out.versions)

    freebayes_vcf = Channel.empty().mix(
                    MERGE_FREEBAYES.out.vcf,
                    bcftools_vcf_out.no_intervals)
                .map{ meta, vcf ->
                    [[
                        id:             meta.id,
                        num_intervals:  meta.num_intervals,
                        variantcaller:  "freebayes"],
                    vcf]
                }
    freebayes_tbi = Channel.empty().mix(
                    MERGE_FREEBAYES.out.tbi,
                    TABIX_VC_FREEBAYES.out.tbi)
                .map{ meta, tbi ->
                    [[
                        id:             meta.id,
                        num_intervals:  meta.num_intervals,
                        variantcaller:  "freebayes"],
                    tbi]
                }

    //Manipulate for payload ingestion
    ch_payload=freebayes_vcf
    .combine(freebayes_tbi)
    .combine(analysis_json)
    .map {metaA,vcf,metaB,tbi,analysis_json ->
    [
        [ id : metaA.id,
            study_id : params.study_id,
            tool : "freebayes"
        ],
        [vcf, tbi], analysis_json]
    }

    //Generate payload
    PAYLOAD_GERMLINEVARIANT(
        ch_payload,
        "",
        "",
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        "freebayes"
    )

    //Gather temporary files
    ch_cleanup=FREEBAYES.out.vcf.map{meta,vcf -> [vcf]}
    .mix(MERGE_FREEBAYES.out.tbi.map{meta,tbi -> [tbi]})
    .mix(BCFTOOLS_SORT.out.vcf.map{meta,vcf -> [vcf]})
    .mix(TABIX_VC_FREEBAYES.out.tbi.map{meta,tbi -> [tbi]})
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

