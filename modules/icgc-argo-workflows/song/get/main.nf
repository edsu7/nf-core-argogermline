process SONG_GET {
    tag "${analysis_id}"
    label 'process_single'

    pod = [secret: workflow.runName + "-secret", mountPath: "/tmp/rdpc_secret"]
    container "${ task.ext.song_container ?: 'ghcr.io/overture-stack/song-client' }:${ task.ext.song_container_version ?: '5.0.2' }"

    if (workflow.containerEngine == "singularity") {
        containerOptions "--bind \$(pwd):/song-client/logs"
    } else if (workflow.containerEngine == "docker") {
        containerOptions "-v \$(pwd):/song-client/logs"
    }

    input:
    tuple val(study_id), val(analysis_id)

    output:
    tuple val(study_id), val(analysis_id), path("*.analysis.json"), emit: analysis_json
    path "*.analysis.json", emit: json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = params.args ?: ''
    def prefix = params.prefix ?: "${analysis_id}"
    def song_url = params.song_url ?: ""
    def accessToken = params.api_token ?: "`cat /tmp/rdpc_secret/secret`"
    def VERSION = params.song_container_version ?: '5.0.2'
    """
    export CLIENT_SERVER_URL=${song_url}
    export CLIENT_STUDY_ID=${study_id}
    export CLIENT_ACCESS_TOKEN=${accessToken}

    sing search -a ${analysis_id} > ${prefix}.analysis.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: ${VERSION}
    END_VERSIONS
    """
}
