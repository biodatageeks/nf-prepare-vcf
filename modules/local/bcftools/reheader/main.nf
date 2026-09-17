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

    # Add only what the original header does not declare already: an ID may be defined once, and
    # CALC_DOSAGE_POLARSBIO passes an input FORMAT/DS definition through rather than replacing it.
    : > header_lines_to_add.txt
    while IFS= read -r line; do
        if [ -z "\$line" ]; then
            continue
        fi
        declared=\$(printf '%s' "\$line" | sed -n 's/^\\(##[A-Za-z]*=<ID=[^,]*,\\).*/\\1/p')
        # index(...) == 1 is a fixed-string match anchored at the start, so an ID quoted inside
        # some other line's Description cannot pass for a declaration.
        if [ -n "\$declared" ] && awk -v d="\$declared" 'index(\$0, d) == 1 { found = 1 } END { exit !found }' header.txt; then
            echo "Original header already declares \${declared}...> - not adding it a second time"
        else
            printf '%s\\n' "\$line" >> header_lines_to_add.txt
        fi
    done <<< '${header_lines_to_add}'

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
