process COMPUTE_GENE_RANGES {

    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/psuszynski/python_tools:1.0.8'

    input:
    tuple val(meta), path(vcf_file), path(script_path)
    val(out_name_part)
    
    output:
    tuple val(meta), path("*_gene_ranges.csv"), emit: gene_ranges
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${script_path} \\
        --input-vcf ${vcf_file} \\
        --output-prefix ${prefix}_${out_name_part} \\
        --threads ${task.cpus}

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
        compute_ranges_py: 0.0.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}_gene_ranges.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS
    """
}