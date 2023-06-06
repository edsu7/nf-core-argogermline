// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { GATK4_HAPLOTYPECALLER                                         } from '../../../modules/nf-core/gatk4/haplotypecaller/main' 
include { GATK4_MERGEVCFS as MERGE_HAPLOTYPECALLER                      } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_CNNSCOREVARIANTS as CNNSCOREVARIANTS                    } from '../../../modules/nf-core/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES as FILTERVARIANTTRANCHES          } from '../../../modules/nf-core/gatk4/filtervarianttranches/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_INDEL                      } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_VCF_SNV                        } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_INDEL_GERMLINEVARIANT  } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { PAYLOAD_GERMLINEVARIANT as PAYLOAD_VCF_SNV_GERMLINEVARIANT    } from '../../../modules/icgc-argo-workflows/payload/germlinevariant/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_INDEL                  } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as SONG_SCORE_UP_VCF_SNV                    } from '../../icgc-argo-workflows/song_score_upload/main'
include { CLEANUP                                                       } from '../../../modules/icgc-argo-workflows/cleanup/main'

workflow GERMLINE_VARIANT_GATK_HAPLOTYPECALLER {
    take:
    cram                            // channel: [mandatory] [meta, cram, crai, interval.bed]
    fasta                           // channel: [mandatory]
    fasta_fai                       // channel: [mandatory]
    dict                            // channel: [mandatory]
    dbsnp                           // channel: []
    dbsnp_tbi
    known_indels
    known_indels_tbi
    known_snps
    known_snps_tbi
    intervals_bed_combined
    analysis_json
    versions

    main:
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(versions)

    //Run Caller
    GATK4_HAPLOTYPECALLER(
           cram,
           fasta,
           fasta_fai,
           dict,
           dbsnp,
           dbsnp_tbi
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    //Split according to VCF and TBI by intervals
    GATK4_HAPLOTYPECALLER.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_vcf_branch}

    GATK4_HAPLOTYPECALLER.out.tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_tbi_branch}

    //Merge VCFs if multiple
    MERGE_HAPLOTYPECALLER(
        haplotypecaller_vcf_branch.intervals
        .map{ meta, vcf ->

            new_meta = [
                        id: meta.id,
                        experimentalStrategy : meta.experimentalStrategy,
                        genomeBuild : meta.genomeBuild,
                        tumourNormalDesignation : meta.tumourNormalDesignation,
                        sampleType : meta.sampleType ,
                        gender : meta.gender,
                        study_id : meta.study_id,
                        num_intervals:  meta.num_intervals
                        ]

                [groupKey(new_meta, new_meta.num_intervals), vcf]
            }.groupTuple(),
        dict.map{ it -> [[id:it[0].baseName], it]})
    ch_versions = ch_versions.mix(MERGE_HAPLOTYPECALLER.out.versions)

    //Rename channels to include merged intervals or single interval
    haplotypecaller_vcf = Channel.empty().mix(
        MERGE_HAPLOTYPECALLER.out.vcf,
        haplotypecaller_vcf_branch.no_intervals)

    haplotypecaller_tbi = Channel.empty().mix(
        MERGE_HAPLOTYPECALLER.out.tbi,
        haplotypecaller_tbi_branch.no_intervals)

    //Manipulate input to combine previously declared channels and intervals_bed_combined
    cnn_in = haplotypecaller_vcf.combine(haplotypecaller_tbi).combine(intervals_bed_combined)
    .map{ metaA, vcf, metaB, tbi, metaC, intervals ->
            new_intervals = intervals.simpleName == "no_intervals" ? [] : intervals
            [metaA, vcf, tbi, [], new_intervals]
        }

    //Score variants
    CNNSCOREVARIANTS(
            cnn_in,
            fasta,
            fasta_fai,
            dict,
            [],
            []
        )
    ch_versions = ch_versions.mix(CNNSCOREVARIANTS.out.versions)

    //Manipulate scored outputs
    cnn_out = CNNSCOREVARIANTS.out.vcf.join(CNNSCOREVARIANTS.out.tbi).combine(intervals_bed_combined)
        .map{   metaA, cnn_vcf,cnn_tbi, metaB, intervals ->
            new_intervals = intervals.simpleName == "no_intervals" ? [] : intervals
            [metaA, cnn_vcf, cnn_tbi, new_intervals]
        }

    //Perform filtering
    FILTERVARIANTTRANCHES(
        cnn_out,
        known_indels.concat(known_snps).flatten().unique().collect(),
        known_indels_tbi.concat(known_snps_tbi).flatten().unique().collect(),
        fasta,
        fasta_fai,
        dict
    )
    ch_versions = ch_versions.mix(FILTERVARIANTTRANCHES.out.versions)

    // Recombine VCF and TBI to make a single payload
    select_vcf=FILTERVARIANTTRANCHES.out.vcf
    .combine(FILTERVARIANTTRANCHES.out.tbi)
    .map{ meta, vcf, metaB, tbi ->
        [[
            id: meta.id,
            experimentalStrategy : meta.experimentalStrategy,
            genomeBuild : meta.genomeBuild,
            tumourNormalDesignation : meta.tumourNormalDesignation,
            sampleType : meta.sampleType ,
            gender : meta.gender,
            study_id : meta.study_id,
            num_intervals:  meta.num_intervals,
            tool : "haplotypecaller"
        ],vcf,tbi,[]]
    }

    //Filter for Indels
    SELECT_VCF_INDEL(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_INDEL.out.versions)


    //Filter for SNVs
    SELECT_VCF_SNV(select_vcf)
    ch_versions = ch_versions.mix(SELECT_VCF_SNV.out.versions)

    //Manipulate for payload ingestion
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
    ch_cleanup=GATK4_HAPLOTYPECALLER.out.vcf.map{meta,vcf -> [vcf]}
    .mix(MERGE_HAPLOTYPECALLER.out.vcf.map{meta,vcf -> [vcf]})
    .mix(CNNSCOREVARIANTS.out.vcf.map{meta,vcf -> [vcf]})
    .mix(FILTERVARIANTTRANCHES.out.vcf.map{meta,vcf -> [vcf]})
    .mix(SELECT_VCF_INDEL.out.vcf.map{meta,vcf -> [vcf]})
    .mix(SELECT_VCF_SNV.out.vcf.map{meta,vcf -> [vcf]})
    .mix(PAYLOAD_VCF_INDEL_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
    .mix(PAYLOAD_VCF_SNV_GERMLINEVARIANT.out.payload_files.map{meta,analysis,files -> [analysis]})
    .collect()

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
    versions       = ch_versions           
}