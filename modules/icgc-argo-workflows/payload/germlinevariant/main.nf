process PAYLOAD_GERMLINEVARIANT {
    tag "$meta.id"
    label 'process_single'


    conda "bioconda::multiqc=1.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:  // input, make update as needed
      tuple val(meta), path(files_to_upload), path(metadata_analysis)
      path pipeline_yml
      val tarball


    output:  // output, make update as needed
      tuple val(meta), path("*.payload.json"), path("out/*"), emit: payload_files
      path "versions.yml", emit: versions

    script:
      // add and initialize variables here as needed
      def arg_pipeline_yml = pipeline_yml.name != 'NO_FILE' ? "-p $pipeline_yml" : ''
      """
      main.py \
        -f ${files_to_upload} \
        -a ${metadata_analysis} \
        -b "${meta.genomeBuild}" \
        -w "DNA Seq Germline Variant Workflow" \
        -r ${workflow.runName} \
        -s "${workflow.sessionId}" \
        -v "${workflow.manifest.version}" \
        -t "${meta.tool}" \
        -l "${tarball}" \
        $arg_pipeline_yml

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$(python --version | sed 's/Python //g')
      END_VERSIONS
      """
  }
