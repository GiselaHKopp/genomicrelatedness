process ANGSD_NGSRELATE {
    tag "${meta.id}"
    label 'process_small'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ngsrelate:2.0--hea85c65_0':
          'biocontainers/ngsrelate:2.0--hea85c65_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: result
    tuple val("${task.process}"), val('ngsrelate'), eval('echo "2.0"'), topic: versions, emit: versions_ngsrelate

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = "ngsrelate.results"

    """
    ngsRelate \\
      -h ${vcf} \\
      -O ${prefix}.${suffix} \\
      ${args} \\
      -p ${task.cpus}
    """
}
