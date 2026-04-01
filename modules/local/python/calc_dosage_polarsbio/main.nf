process CALC_DOSAGE_POLARSBIO {

    tag "$meta.id"
    label 'process_1'

    container 'docker.io/psuszynski/python_tools:1.0.11'

    input:
    tuple val(meta), path(vcf_file), path(python_script)
    val(calc_ds_min_gq)
    val(out_name_part)
    
    output:
    tuple val(meta), path("*_${out_name_part}.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    python3 ${python_script} \\
        --input-vcf-path ${vcf_file} \\
        --output-vcf-path ${prefix}_${out_name_part}.vcf.bgz \\
        --calc-ds-min-gq ${calc_ds_min_gq}
    
    mv ${prefix}_${out_name_part}.vcf.bgz ${prefix}_${out_name_part}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polars-bio: 0.26.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polars-bio: 0.26.1
    END_VERSIONS
    """
}