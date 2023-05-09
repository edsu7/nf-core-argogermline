// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { CNVKIT_BATCH                             } from '../../../modules/nf-core/cnvkit/batch/main'

workflow GERMLINE_VARIANT_CNVKIT {

    take:
        cram
        fasta
        fasta_fai
        target
        reference
        panel_of_normals
        analysis_json
        versions
    main:

        ch_versions = Channel.empty()

        CNVKIT_BATCH(cram, fasta, fasta_fai, target, reference,panel_of_normals)

        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

    emit:
    bed = CNVKIT_BATCH.out.bed
    cnn = CNVKIT_BATCH.out.cnn
    cnr = CNVKIT_BATCH.out.cnr
    cns = CNVKIT_BATCH.out.cns
    pdf = CNVKIT_BATCH.out.pdf
    png = CNVKIT_BATCH.out.png
    versions = ch_versions                     // channel: [ versions.yml ]
}

