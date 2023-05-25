import groovy.json.JsonSlurper

process GET_META {
    //tag "$meta.id"
    label 'process_single'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jq:1.6' : 
        'ghcr.io/overture-stack/song-client' }"

    input:
    path(analysis_json)

    output:
    //val env(gender), emit: gender
    //val env(tumourNormalDesignation) , emit: designation
    //val env(sampleType) , emit: sampleType
    //val env(genomeBuild), emit : genomeBuild
    val experimentalStrategy
    when:
    task.ext.when == null || task.ext.when

    script:
    experimentalStrategy = new JsonSlurper().parseText(analysis_json.text).get('experiment').get('experimental_strategy')
    // """
    // gender=\$(cat ${analysis_json} | jq '.samples[0].donor.gender' | sed 's/"//g' )
    // tumourNormalDesignation=\$(cat ${analysis_json} | jq '.samples[0].sampleType' | sed 's/"//g' )
    // sampleType=\$(cat ${analysis_json} | jq '.samples[0].specimen.tumourNormalDesignation' | sed 's/"//g')
    // genomeBuild=\$(cat ${analysis_json} | jq '.workflow.genome_build' | sed 's/"//g' )
    // experimentalStrategy=\$(cat ${analysis_json} | jq '.experiment.experimental_strategy' |  sed 's/"//g' )

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     local: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    // END_VERSIONS
    // """
}
