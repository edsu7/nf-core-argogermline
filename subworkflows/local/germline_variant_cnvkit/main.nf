// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { CNVKIT_BATCH                                     } from '../../../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_EXPORT                                    } from '../../../modules/nf-core/cnvkit/export/main'
include { TABIX_BGZIPTABIX                                 } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VARIANTS      } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_SUPPLEMENT    } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD as UP_VARIANTS                 } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UP_SUPPLEMENT               } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                                          } from '../../../modules/icgc-argo-workflows/cleanup/main'

workflow GERMLINE_VARIANT_CNVKIT {

    take:
        cram
        fasta
        fasta_fai
        target
        analysis_json
        versions
    main:

        ch_versions = Channel.empty()

        //Run Caller
        CNVKIT_BATCH(cram, fasta, fasta_fai, target, [],false)
        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

        //Three possible files are generate : 
        //".sorted.bintest.cns",".sorted.call.cns",".sorted.cns"
        //Use ".sorted.call.cns"
        ch_export=CNVKIT_BATCH.out.cns.map{ meta,files -> [meta,files.sort().get(1)]}

        //Convert cns to vcf
        CNVKIT_EXPORT(ch_export)
        ch_versions = ch_versions.mix(CNVKIT_EXPORT.out.versions)

        //Compress and index
        TABIX_BGZIPTABIX(CNVKIT_EXPORT.out.output)
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)


        //Manipulate for payload ingestion
        ch_supplement_payload=CNVKIT_BATCH.out.cns\
        .combine(CNVKIT_BATCH.out.cnr.map{meta,files->[files]})
        .combine(CNVKIT_BATCH.out.cnn.map{meta,files->[files]})
        .combine(CNVKIT_BATCH.out.bed.map{meta,files->[files]})
        .combine(CNVKIT_BATCH.out.png.map{meta,files->[files]})
        .combine(CNVKIT_BATCH.out.pdf.map{meta,files->[files]})
        .combine(analysis_json)
        .map{
            meta,cns,cnr,cnn,bed,png,pdf,analysis_json ->
            [
                [
                id : meta.id,
                experimentalStrategy : meta.experimentalStrategy,
                genomeBuild : meta.genomeBuild,
                tumourNormalDesignation : meta.tumourNormalDesignation,
                sampleType : meta.sampleType ,
                gender : meta.gender,
                study_id : params.study_id,
                dataType : "CNV",
                tool : "cnvkit"
                ]
                ,[cns,cnr,cnn,bed,png,pdf].flatten().collect(),analysis_json
            ]
        }

        ch_payload=TABIX_BGZIPTABIX.out.gz_tbi.combine(analysis_json)
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
                    dataType : "CNV",
                    tool : "cnvkit"
                    ]
                    ,[vcf, tbi],analysis_json]
                }

        //Generate payload
        PAYLOAD_VARIANTS(
            ch_payload,
            ch_versions.unique().collectFile(name: 'collated_versions.yml'),
            false
        )

        PAYLOAD_SUPPLEMENT(
            ch_supplement_payload,
            ch_versions.unique().collectFile(name: 'collated_versions.yml'),
            true
        )

        //Gather temporary files
        ch_cleanup=ch_export.map{meta,cns -> [cns]}
        .mix(CNVKIT_EXPORT.out.output.map{meta,output -> [output]})
        .mix(TABIX_BGZIPTABIX.out.gz_tbi.map { meta,vcf,tbi -> [vcf]})
        .mix(PAYLOAD_VARIANTS.out.payload_files.map{meta,analysis,files -> [analysis]})
        .mix(PAYLOAD_SUPPLEMENT.out.payload_files.map{meta,analysis,files -> [analysis]})
        .collect()

        //If Local is true, will be published into "output_dir" directory
        if (params.local==false){
            //Upload variants
            //UP_VARIANT(PAYLOAD_VARIANTS.out.payload_files)
            //UP_SUPPLEMENT(PAYLOAD_SUPPLEMENT.out.payload_files)
            if (params.cleanup){
                CLEANUP(
                    ch_cleanup.collect(),
                    UP_VARIANT.out.analysis_id.mix(UP_SUPPLEMENT.out.analysis_id).collect()
                    )
            }
        } else {
        if (params.cleanup){
            CLEANUP(
                ch_cleanup.collect(),
                PAYLOAD_VARIANTS.out.payload_files.map{meta,analysis,files -> [analysis]}
                .mix(PAYLOAD_SUPPLEMENT.out.payload_files.map{meta,analysis,files -> [analysis]}).collect()
                )
            }       
        }



    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}

