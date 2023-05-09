// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { TIDDIT_SV                         } from '../../../modules/nf-core/tiddit/sv/main'
include { TABIX_BGZIPTABIX                  } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { PAYLOAD_GERMLINEVARIANT           } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD                 } from '../../icgc-argo-workflows/song_score_upload/main'
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
                [ id : meta.id,
                  study_id : params.study_id,
                  tool : "tiddit"
                ]
                ,[vcf, tbi],analysis_json]
            }.set{ch_payload}

    //Generate payload
    PAYLOAD_GERMLINEVARIANT(
        ch_payload,
        "",
        "",
        ch_versions.unique().collectFile(name: 'collated_versions.yml'),
        "tiddit"
    )

    //Gather temporary files
    ch_cleanup=TIDDIT_SV.out.vcf.combine(TABIX_BGZIPTABIX.out.gz_tbi)
        .map{ metaA,vcf,metaB,vcf_gz,vcf_gz_tbi ->
            [vcf,vcf_gz,vcf_gz_tbi]
        }
       .collect()
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

