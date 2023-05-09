// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { GATK4_GENOMICSDBIMPORT                                } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { BCFTOOLS_SORT                                         } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX as TABIX                                  } from '../../../modules/nf-core/tabix/tabix/main'
include { GATK4_GENOTYPEGVCFS                                   } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS                                       } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_INDEL} from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_SNP  } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDEL              } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNP                } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_MERGEVCFS as MERGE_GENOTYPEGVCFS                } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_VQSR                         } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow GERMLINE_VARIANT_ANNOTATION {

    take:
    cram                            // channel: [mandatory] [meta, cram, crai, interval.bed]
    fasta                           // channel: [mandatory]
    fasta_fai                       // channel: [mandatory]
    dict                            // channel: [mandatory]
    dbsnp                           // channel: []
    dbsnp_tbi
    known_sites_indels
    known_sites_indels_tbi
    known_sites_snps
    known_sites_snps_tbi

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    // // Figure out if using intervals or not
    GATK4_HAPLOTYPECALLER.out.vcf.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_vcf_branch}

    GATK4_HAPLOTYPECALLER.out.tbi.branch{
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }.set{haplotypecaller_tbi_branch}

    // // Merge VCFs
    genotype_gvcf_to_call = Channel.empty().mix(GATK4_HAPLOTYPECALLER.out.vcf
                                                .join(GATK4_HAPLOTYPECALLER.out.tbi)
                                                .join(cram).map{ meta, vcf, tbi, cram, crai, intervals, dragstr_model ->
                                                        [ meta, vcf, tbi, intervals ]
                                                })

    // // Prepare input for GenomicsDBImport
    gendb_input = genotype_gvcf_to_call.map{
        meta, gvcf, tbi, intervals->
            new_meta = [
                        id:             "joint_variant_calling",
                        intervals_name: meta.intervals_name,
                        num_intervals:  meta.num_intervals
                    ]

            [ new_meta, gvcf, tbi, intervals ]

        }.groupTuple(by:[0,3]).map{ new_meta, gvcf, tbi, intervals ->
            [ new_meta, gvcf, tbi, intervals, [], [] ]
        }

    // //Convert VCFs into genomicsDB
    GATK4_GENOMICSDBIMPORT ( gendb_input, false, false, false )

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb.map{
        meta, genomicsdb ->
            [meta, genomicsdb, [], [], []]
        }

    // //Joint Genotyping via GenoTypeGVCFs
    GATK4_GENOTYPEGVCFS(genotype_input, fasta, fasta_fai, dict, dbsnp, dbsnp_tbi)

    // //Sort each VCFs
    BCFTOOLS_SORT(GATK4_GENOTYPEGVCFS.out.vcf)

    // // Split sorted VCFs if intervals
    vcfs_sorted_input = BCFTOOLS_SORT.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // // Gather into tuple
    vcfs_sorted_input_no_intervals =  vcfs_sorted_input.no_intervals.map{ meta , vcf ->

                            [[  id:             "joint_variant_calling",
                                num_intervals:  meta.num_intervals,
                                patient:        "all_samples",
                                variantcaller:  "haplotypecaller"
                            ] , vcf ]
    }

    // //Index VCF files
    TABIX(vcfs_sorted_input_no_intervals)

    // //Merged scattered VCFs and gather into tuple
    MERGE_GENOTYPEGVCFS(vcfs_sorted_input.intervals.map{meta, vcf ->
                            [
                                    [
                                        id: "joint_variant_calling",
                                        num_intervals: meta.num_intervals,
                                        patient: "all_samples",
                                        variantcaller: "haplotypecaller",
                                    ],
                            vcf]
                        }.groupTuple()
                        ,dict)

    // //Input channel for annotation
    vqsr_input = Channel.empty().mix(
        MERGE_GENOTYPEGVCFS.out.vcf.join(MERGE_GENOTYPEGVCFS.out.tbi),
        vcfs_sorted_input_no_intervals.join(TABIX.out.tbi)
    )

    // //Group resources for SNP and INDEL respectively
    snp_resource_labels   = Channel.empty().mix(known_snps_vqsr,dbsnp_vqsr).collect()
    indel_resource_labels = Channel.empty().mix(known_indels_vqsr,dbsnp_vqsr).collect()
    
    // //Recalibrate SNP
    VARIANTRECALIBRATOR_SNP(
        vqsr_input,
        dbsnp,
        dbsnp_tbi,
        snp_resource_labels,
        fasta,
        fasta_fai,
        dict)
    // //Recalibrate INDEL
    VARIANTRECALIBRATOR_INDEL(
        vqsr_input,
        dbsnp,
        dbsnp_tbi,
        indel_resource_labels,
        fasta,
        fasta_fai,
        dict)

    // // Gather recalibrated SNPs into tuple
    vqsr_input_snp   = vqsr_input.join( VARIANTRECALIBRATOR_SNP.out.recal)
                                .join( VARIANTRECALIBRATOR_SNP.out.idx)
                                .join( VARIANTRECALIBRATOR_SNP.out.tranches)
                                .map{ meta, vcf, tbi, recal, index, tranche ->

                                            new_meta = [
                                                        id:             "recalibrated_joint_variant_calling",
                                                        num_intervals:  meta.num_intervals,
                                                        patient:        "all_samples",
                                                        variantcaller:  "haplotypecaller",
                                                    ]

                                            [new_meta, vcf, tbi, recal, index, tranche]
                                        }
    // // Gather recalibrated INDELs into tuple
    vqsr_input_indel = vqsr_input.join( VARIANTRECALIBRATOR_INDEL.out.recal).join(
                                        VARIANTRECALIBRATOR_INDEL.out.idx).join(
                                        VARIANTRECALIBRATOR_INDEL.out.tranches).map{ meta, vcf, tbi, recal, index, tranche ->

                                                new_meta = [
                                                            id:             "recalibrated_joint_variant_calling",
                                                            num_intervals:  meta.num_intervals,
                                                            patient:        "all_samples",
                                                            variantcaller:  "haplotypecaller"
                                                ]

                                            [new_meta, vcf, tbi, recal, index, tranche]
                                        }

    // //Apply VSQR filtering for SNPs and INDEL
    GATK4_APPLYVQSR_SNP(vqsr_input_snp,
                        fasta,
                        fasta_fai,
                        dict )

    GATK4_APPLYVQSR_INDEL(vqsr_input_indel,
                        fasta,
                        fasta_fai,
                        dict )
    
    // //Merge VSQR filtered VCFs

    vqsr_snp_vcf = GATK4_APPLYVQSR_SNP.out.vcf
    vqsr_indel_vcf = GATK4_APPLYVQSR_INDEL.out.vcf

    MERGE_VQSR(
        vqsr_snp_vcf.mix(vqsr_indel_vcf).groupTuple(),
        dict
    )

    ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
    ch_versions = ch_versions.mix(VARIANTRECALIBRATOR_SNP.out.versions)
    ch_versions = ch_versions.mix(GATK4_APPLYVQSR_SNP.out.versions)

    emit:
    versions       = ch_versions   
    //vcf   = Channel.empty().mix( vcfs_sorted_input_no_intervals, MERGE_GENOTYPEGVCFS.out.vcf, MERGE_VQSR.out.vcf)  // channel: [ val(meta), [ vcf ] ]
    //index = Channel.empty().mix( TABIX.out.tbi, MERGE_GENOTYPEGVCFS.out.tbi, MERGE_VQSR.out.tbi) 
}

