process BCFTOOLS_REHEADER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(original_vcf), path(vcf_with_incorrect_header)
    val(header_lines_to_add)
    val(out_name_part)

    output:
    tuple val(meta), path("*_${out_name_part}.vcf.gz"),     emit: vcf
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("${original_vcf}" == "${prefix}_${out_name_part}.vcf.gz" || "${vcf_with_incorrect_header}" == "${prefix}_${out_name_part}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    bcftools view -h ${original_vcf} > header.txt

    head -n -1 header.txt > header_without_last_line.txt
    tail -n 1 header.txt > header_last_line.txt

    # remove last line of header_without_last_line.txt if it starts with '##bcftools_viewCommand=view -h chr22_filtered.vcf.gz;' because this lines contans a date and this makes the modules/local/bcftools/reheader/tests/main.nf.test snapshot martch fail
    if tail -n 1 header_without_last_line.txt | grep -q '^##bcftools_viewCommand=view -h'; then
        head -n -1 header_without_last_line.txt > header_without_last_line_tmp.txt
        mv header_without_last_line_tmp.txt header_without_last_line.txt
    fi

    echo '${header_lines_to_add}' > header_lines_to_add.txt
    cat header_without_last_line.txt header_lines_to_add.txt header_last_line.txt > header_final.txt

    bcftools reheader \\
        ${args} \\
        -h header_final.txt \\
        --output ${prefix}_${out_name_part}.vcf.gz \\
        --threads ${task.cpus} \\
        ${vcf_with_incorrect_header}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${out_name_part}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
