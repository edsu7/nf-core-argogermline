// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { CREATE_INTERVALS_BED                                          } from '../../../modules/local/create_intervals_bed/main'
include { TABIX_BGZIPTABIX                                              } from '../../../modules/nf-core/tabix/bgziptabix/main'
workflow PREPARE_INTERVALS {

    take:
        intervals_bed

    main:
        ch_versions = Channel.empty()

        ch_intervals                     = Channel.empty() // List of bed files, one for each region
        ch_intervals_bed_gz_tbi          = Channel.empty() // List of bed.gz, bed,gz.tbi, one for each region
        ch_intervals_combined            = Channel.empty() // Bed file containing all intervals

        ch_intervals_combined = intervals_bed.map{it -> [[id:it.baseName], it] }
        
        ch_intervals = CREATE_INTERVALS_BED(intervals_bed.collect()).bed

        ch_intervals = ch_intervals.flatten()
            .map{ intervalFile ->
                def duration = 0.0
                for (line in intervalFile.readLines()) {
                    final fields = line.split('\t')
                    if (fields.size() >= 5) duration += fields[4].toFloat()
                    else {
                        start = fields[1].toInteger()
                        end = fields[2].toInteger()
                        duration += (end - start) / params.nucleotides_per_second
                    }
                }
                [duration, intervalFile]
            }.toSortedList({ a, b -> b[0] <=> a[0] })
            .flatten().collate(2)
            .map{duration, intervalFile -> intervalFile}
            .collect().map{ it ->
                    [it, it.size() ] // Adding number of intervals as elements
                }.transpose()


        tabix_in = ch_intervals.map{ file, num_intervals -> [[id:file.baseName], file] }

        TABIX_BGZIPTABIX(tabix_in)
        ch_intervals_bed_gz_tbi = TABIX_BGZIPTABIX.out.gz_tbi.map{ meta, bed, tbi -> [bed, tbi ]}.toList().map{
                                        it ->
                                        [it, it.size()] // Adding number of intervals as elements
                                    }.transpose()
    
        ch_versions = ch_versions.mix(CREATE_INTERVALS_BED.out.versions)
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:

        intervals_bed               = ch_intervals                                           // path: intervals.bed, num_intervals                        [intervals split for parallel execution]
        intervals_bed_gz_tbi        = ch_intervals_bed_gz_tbi                                // path: target.bed.gz, target.bed.gz.tbi, num_intervals     [intervals split for parallel execution]
        intervals_bed_combined      = ch_intervals_combined                                  // path: intervals.bed                        [all intervals in one file]
        versions                    = ch_versions                                            // channel: [ versions.yml ]

}

