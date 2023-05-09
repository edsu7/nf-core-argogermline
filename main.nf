#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/argogermline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/argogermline

    Website: https://nf-co.re/argogermline
    Slack  : https://nfcore.slack.com/channels/argogermline
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.study_id             = WorkflowMain.getGenomeAttribute(params, 'study_id')
params.analysis_id          = WorkflowMain.getGenomeAttribute(params, 'analysis_id')

params.api_token            = WorkflowMain.getGenomeAttribute(params, 'api_token')
params.score_url_upload     = WorkflowMain.getGenomeAttribute(params, 'score_url_upload')
params.song_url_upload      = WorkflowMain.getGenomeAttribute(params, 'song_url_upload')
params.score_url_download   = WorkflowMain.getGenomeAttribute(params, 'score_url_download')
params.song_url_download    = WorkflowMain.getGenomeAttribute(params, 'song_url_download')
params.score_url            = WorkflowMain.getGenomeAttribute(params, 'score_url')
params.song_url             = WorkflowMain.getGenomeAttribute(params, 'song_url')

params.dict                 = WorkflowMain.getGenomeAttribute(params, 'dict')
params.fasta_fai            = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.fasta                = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.known_snps           = WorkflowMain.getGenomeAttribute(params, 'known_snps')
params.known_snps_tbi       = WorkflowMain.getGenomeAttribute(params, 'known_snps_tbi')
params.dbsnp                = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi            = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.bwa                  = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.intervals_bed        = WorkflowMain.getGenomeAttribute(params, 'intervals_bed')
params.intervals_bed_gz     = WorkflowMain.getGenomeAttribute(params, 'intervals_bed_gz')
params.intervals_bed_gz_tbi = WorkflowMain.getGenomeAttribute(params, 'intervals_bed_gz_tbi')
params.known_indels         = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi     = WorkflowMain.getGenomeAttribute(params, 'known_indels_tbi')

params.tools                = WorkflowMain.getGenomeAttribute(params, 'tools')
params.outdir               = WorkflowMain.getGenomeAttribute(params, 'outdir')
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ARGOGERMLINE } from './workflows/argogermline'

//
// WORKFLOW: Run main nf-core/argogermline analysis pipeline
//
workflow NFCORE_ARGOGERMLINE {
    ARGOGERMLINE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_ARGOGERMLINE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
