include { SONG_SCORE_DOWNLOAD as NORMAL_SONG_SCORE_DOWNLOAD             } from '../../icgc-argo-workflows/song_score_download/main'

include { CRAM_RECALIBRATE as NORMAL_GATK4_RECALIBRATE                  } from '../cram_recalibrate/main'

include { GERMLINE_VARIANT_STRELKA                                      } from '../germline_variant_strelka/main'
include { GERMLINE_VARIANT_DEEPVARIANT                                  } from '../germline_variant_deepvariant/main'
include { GERMLINE_VARIANT_TIDDIT                                       } from '../germline_variant_tiddit/main'
include { GERMLINE_VARIANT_MANTA                                        } from '../germline_variant_manta/main'
include { GERMLINE_VARIANT_CNVKIT                                       } from '../germline_variant_cnvkit/main'
include { GERMLINE_VARIANT_FREEBAYES                                    } from '../germline_variant_freebayes/main'
include { GERMLINE_VARIANT_GATK_HAPLOTYPECALLER                         } from '../germline_variant_gatk_haplotypecaller/main'
include { PREPARE_INTERVALS                                             } from '../prepare_intervals/main'
include { BUILD_INTERVALS                                               } from '../../../modules/local/build_intervals/main'
include { CLEANUP                                                       } from '../../../modules/icgc-argo-workflows/cleanup/main'

fasta                = params.fasta                ? Channel.fromPath(params.fasta).collect()                              : Channel.empty()
fasta_fai            = params.fasta_fai            ? Channel.fromPath(params.fasta_fai).collect()                          : Channel.empty()
germline_resource    = params.germline_resource    ? Channel.fromPath(params.germline_resource).collect()                  : Channel.value([]) //Mutec2 does not require a germline resource, so set to optional input
known_indels         = params.known_indels         ? Channel.fromPath(params.known_indels).collect()                       : Channel.empty()
known_indels_tbi     = params.known_indels_tbi     ? Channel.fromPath(params.known_indels_tbi).collect()                   : Channel.empty()
known_snps           = params.known_snps           ? Channel.fromPath(params.known_snps).collect()                         : Channel.empty()
known_snps_tbi       = params.known_snps_tbi       ? Channel.fromPath(params.known_snps_tbi).collect()                     : Channel.empty()
dbsnp                = params.dbsnp                ? Channel.fromPath(params.dbsnp).collect()                              : Channel.empty()
dbsnp_tbi            = params.dbsnp_tbi            ? Channel.fromPath(params.dbsnp_tbi).collect()                          : Channel.empty()
dragstr_model        = params.dragstr_model        ? Channel.fromPath(params.dragstr_model).collect()                      : Channel.empty()
dict                 = params.dict                 ? Channel.fromPath(params.dict).collect()                               : Channel.empty()
bwa                  = params.bwa                  ? Channel.fromPath(params.bwa).collect()                                : Channel.value([])
wes                  = params.wes                  ? params.wes                                                            : false

workflow GERMLINE_VARIANTS {

    take:
    // TODO nf-core: edit input (take) channels
        study_id
        analysis_id
        tools

    main:

    ch_versions = Channel.empty()
    cumulative_versions = Channel.empty()

    //ITEMS LEFT TO DO
    //Recalibrate supplemental files
    //Tiddit supplemental files
    //CNVKit convert to VCF?
    //CNVKIt supplemental files

    //Exome conf
    //manta,strelka,deepvariant,cnvkit

    //Gender conf
    //cnvkit

    //Download Files
    NORMAL_SONG_SCORE_DOWNLOAD(tuple([study_id,analysis_id]))
    ch_versions=ch_versions.mix(NORMAL_SONG_SCORE_DOWNLOAD.out.versions)
    cumulative_versions=cumulative_versions.mix(NORMAL_SONG_SCORE_DOWNLOAD.out.versions)

    NORMAL_SONG_SCORE_DOWNLOAD.out.analysis_json.map{
        it -> 
        [
            experimentalStrategy : new groovy.json.JsonSlurper().parse(it).get('experiment').get('experimental_strategy'),
            genomeBuild : new groovy.json.JsonSlurper().parse(it).get('workflow').get('genome_build'),
            tumourNormalDesignation : new groovy.json.JsonSlurper().parse(it).get('samples').get(0).get("specimen").get("tumourNormalDesignation"),
            sampleType : new groovy.json.JsonSlurper().parse(it).get("samples").get(0).get("sampleType"),
            gender : new groovy.json.JsonSlurper().parse(it).get("samples").get(0).get("donor").get("gender"),
            id : new groovy.json.JsonSlurper().parse(it).get("analysisId")
            ]
        }.set{ ch_meta }
    
    //Generate intervals
    if (params.no_intervals){
        intervals_bed = Channel.fromPath(params.intervals_bed, checkIfExists: true).map{it -> [it,1]}

        intervals_bed_gz_tbi=Channel.fromPath(params.intervals_bed_gz, checkIfExists: true)
        .combine(Channel.fromPath(params.intervals_bed_gz_tbi, checkIfExists: true))
        .map{bed_gz,bed_gz_tbi -> [[bed_gz,bed_gz_tbi],1]
        }

        intervals_bed_combined = Channel.fromPath(params.intervals_bed, checkIfExists: true)
        .map{it -> [[id:it.SimpleName],it]}
    } else if (params.intervals_bed){
        PREPARE_INTERVALS(Channel.fromPath(params.intervals_bed,checkIfExists: true).collect())
        intervals_bed = PREPARE_INTERVALS.out.intervals_bed
        intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi
        intervals_bed_combined = PREPARE_INTERVALS.out.intervals_bed_combined
    } else {
        BUILD_INTERVALS(fasta_fai.map{it->[[id:it.SimpleName],it]})
        PREPARE_INTERVALS(BUILD_INTERVALS.out.bed)
        intervals_bed = PREPARE_INTERVALS.out.intervals_bed
        intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi
        intervals_bed_combined = PREPARE_INTERVALS.out.intervals_bed_combined
    }

    //Recalibrate
    if (params.tools.split(',').contains('recalibrate')){
        ch_normal_cram_for_recalibration=NORMAL_SONG_SCORE_DOWNLOAD.out.files.combine(ch_meta)
            .map{cram,crai,meta -> 
            [meta,cram,crai]
            }

        NORMAL_GATK4_RECALIBRATE(
            ch_normal_cram_for_recalibration, //tuple val(meta), path(input), path(input_index), path (target_bed), path (target_bed_tbi)
            dict,
            fasta,
            fasta_fai,
            Channel.value([file("NO_FILE"),0]),
            known_snps,
            known_snps_tbi
        )

        ch_cram_variant_calling_normal = NORMAL_GATK4_RECALIBRATE.out.cram
            .map { meta,cram,crai->
                [meta, cram, crai]
        }
        ch_versions=ch_versions.mix(NORMAL_GATK4_RECALIBRATE.out.versions)
        cumulative_versions=cumulative_versions.mix(NORMAL_GATK4_RECALIBRATE.out.versions)
    } else {
        ch_cram_variant_calling_normal = NORMAL_SONG_SCORE_DOWNLOAD.out.files.combine(ch_meta)
            .map { cram,crai,meta->
                [meta, cram, crai]
            }
    }

    //DEEP VARIANT
    if (params.tools.split(',').contains('deepvariant')){
        ch_normal_deepvariant=ch_cram_variant_calling_normal.combine(intervals_bed)//.combine(Channel.value([file("NO_FILE"),0]))
            .map{meta,cram,crai,intervals,num_intervals ->
                    intervals_new = num_intervals == 0 ? [] : intervals
                    [
                        [
                            id : meta.id,
                            experimentalStrategy : meta.experimentalStrategy,
                            genomeBuild : meta.genomeBuild,
                            tumourNormalDesignation : meta.tumourNormalDesignation,
                            sampleType : meta.sampleType ,
                            gender : meta.gender,
                            num_intervals : num_intervals
                ],  cram, crai, intervals_new]
            }

        GERMLINE_VARIANT_DEEPVARIANT(
            ch_normal_deepvariant,
            dict,
            fasta.map{ it -> [[id:it[0].baseName], it] },
            fasta_fai.map{ it -> [[id:it[0].baseName], it] },
            NORMAL_SONG_SCORE_DOWNLOAD.out.analysis_json,
            ch_versions
        )
        cumulative_versions=cumulative_versions.mix(GERMLINE_VARIANT_DEEPVARIANT.out.versions)
    }

    //CNVKIT
    if (params.tools.split(',').contains('cnvkit')){
        ch_normal_cnvkit = ch_cram_variant_calling_normal
            .map{meta,cram,crai->
                [meta,[],cram]
        }

        GERMLINE_VARIANT_CNVKIT(
            ch_normal_cnvkit, 
            fasta,
            fasta_fai,
            intervals_bed_combined.map{ meta,intervals -> intervals },
            NORMAL_SONG_SCORE_DOWNLOAD.out.analysis_json,
            ch_versions
        )
        cumulative_versions=cumulative_versions.mix(GERMLINE_VARIANT_CNVKIT.out.versions)
    }

    //TIDDIT N
    if (params.tools.split(',').contains('tiddit')){
        ch_normal_tiddit=ch_cram_variant_calling_normal
            .map{meta,cram,crai->
                [meta,cram,crai]
            }

        GERMLINE_VARIANT_TIDDIT(
            ch_normal_tiddit,
            fasta.map{ it -> [[id:it[0].baseName], it] },
            bwa.map{ it -> [[id:it[0].baseName], it] },
            NORMAL_SONG_SCORE_DOWNLOAD.out.analysis_json,
            ch_versions
            )

        cumulative_versions=cumulative_versions.mix(GERMLINE_VARIANT_TIDDIT.out.versions)
    }
 
    // MANTA N
    if (params.tools.split(',').contains('manta')){
        ch_normal_manta=ch_cram_variant_calling_normal.combine(intervals_bed_gz_tbi)
            .map{meta, cram, crai, intervals_bed_gz_tbi, num_intervals ->
                [
                    [
                        id : meta.id,
                        experimentalStrategy : meta.experimentalStrategy,
                        genomeBuild : meta.genomeBuild,
                        tumourNormalDesignation : meta.tumourNormalDesignation,
                        sampleType : meta.sampleType ,
                        gender : meta.gender,
                        num_intervals : num_intervals
                ], cram, crai, intervals_bed_gz_tbi[0],intervals_bed_gz_tbi[1]]
            }

        GERMLINE_VARIANT_MANTA(
            ch_normal_manta,
            dict,
            fasta,
            fasta_fai,
            NORMAL_SONG_SCORE_DOWNLOAD.out.analysis_json,
            ch_versions
            )
        cumulative_versions=cumulative_versions.mix(GERMLINE_VARIANT_MANTA.out.versions)
    }
    // STRELKA N
    if (params.tools.split(',').contains('strelka')){
        ch_normal_strelka=ch_cram_variant_calling_normal.combine(intervals_bed_gz_tbi)
        .map{meta, cram, crai, intervals_bed_gz_tbi, num_intervals ->
                [
                    [
                        id : meta.id,
                        experimentalStrategy : meta.experimentalStrategy,
                        genomeBuild : meta.genomeBuild,
                        tumourNormalDesignation : meta.tumourNormalDesignation,
                        sampleType : meta.sampleType ,
                        gender : meta.gender,
                        num_intervals : num_intervals
                ], cram, crai, intervals_bed_gz_tbi[0], intervals_bed_gz_tbi[1]]
        }

        GERMLINE_VARIANT_STRELKA(
        ch_normal_strelka,
        dict,
        fasta,
        fasta_fai,
        NORMAL_SONG_SCORE_DOWNLOAD.out.analysis_json,
        ch_versions
        )

        //ch_vcf=ch_vcf.mix(GERMLINE_VARIANT_STRELKA.out.vcf)
        cumulative_versions=cumulative_versions.mix(GERMLINE_VARIANT_STRELKA.out.versions)
    }


    //GATK_HAPLOTYPECALLER N
    if (params.tools.split(',').contains('haplotypecaller')){
        ch_normal_gatk_halpotypecaller = ch_cram_variant_calling_normal.combine(intervals_bed)
            .map{ meta, cram, crai, intervals_bed,num_intervals ->
                [
                    [
                        id : meta.id,
                        experimentalStrategy : meta.experimentalStrategy,
                        genomeBuild : meta.genomeBuild,
                        tumourNormalDesignation : meta.tumourNormalDesignation,
                        sampleType : meta.sampleType ,
                        gender : meta.gender,
                        num_intervals : num_intervals
                ], 
                cram, crai, intervals_bed, []]
        }
    
        GERMLINE_VARIANT_GATK_HAPLOTYPECALLER(
            ch_normal_gatk_halpotypecaller,
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            known_indels,
            known_indels_tbi,
            known_snps,
            known_snps_tbi,  
            intervals_bed_combined,
            NORMAL_SONG_SCORE_DOWNLOAD.out.analysis_json,
            ch_versions
        )

        //ch_vcf=ch_vcf.mix(GERMLINE_VARIANT_GATK_HAPLOTYPECALLER.out.vcf)
        cumulative_versions=cumulative_versions.mix(GERMLINE_VARIANT_GATK_HAPLOTYPECALLER.out.versions)
    }

    //FREEBAYES N
    if (params.tools.split(',').contains('freebayes')){
        ch_normal_freebayes = ch_cram_variant_calling_normal.combine(intervals_bed)
            .map{ meta, cram, crai, intervals_bed,num_intervals ->
                [
                    [
                        id : meta.id,
                        experimentalStrategy : meta.experimentalStrategy,
                        genomeBuild : meta.genomeBuild,
                        tumourNormalDesignation : meta.tumourNormalDesignation,
                        sampleType : meta.sampleType ,
                        gender : meta.gender,
                        num_intervals : num_intervals
                ], 
                cram, crai, [], [], intervals_bed]
        }
    
        GERMLINE_VARIANT_FREEBAYES(
            ch_normal_freebayes,
            dict,
            fasta,
            fasta_fai,
            NORMAL_SONG_SCORE_DOWNLOAD.out.analysis_json,
            ch_versions
        )

        cumulative_versions=cumulative_versions.mix(GERMLINE_VARIANT_FREEBAYES.out.versions)
    }

    // Recal CRAM clean up
    if (params.cleanup){
        CLEANUP(ch_cram_variant_calling_normal.map{meta,cram,crai -> [cram,crai]}.collect(),cumulative_versions.collect())
    }  

    emit:
        versions = cumulative_versions                     // channel: [ versions.yml ]
}

