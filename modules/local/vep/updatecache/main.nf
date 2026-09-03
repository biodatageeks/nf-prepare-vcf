process VEP_UPDATECACHE {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:release_113.4':
        'docker.io/psuszynski/ensembl-vep:113.4.3' }"

    input:
    tuple val(meta), path(vep_cache)
    val(species)
    val(input_args)
    val(vep_cache_url)
    tuple val(vep_ref_fasta_url), val(vep_fasta_path)
    // input_ref_fasta / input_ref_fasta_index are optional ([] when not supplied): when given,
    // that FASTA is used as the reference as-is and only a missing .fai is built here.
    tuple path(input_ref_fasta), path(input_ref_fasta_index), val(default_ref_fasta_url)

    output:
    path("${vep_cache}/${species}"), emit: cachesubdir
    path("ref_fasta.fa"), emit: ref_fasta
    path("ref_fasta.fa.fai"), emit: ref_fasta_index
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ ! -d "${vep_cache}/$species" ]; then
        if [ -z "${vep_cache_url}" ]; then
            echo "VEP cache URL is empty, running INSTALL.pl to download the cache for species ${species}"
            perl /opt/vep/src/ensembl-vep/INSTALL.pl --CACHEDIR ${vep_cache} --SPECIES $species $args $input_args
            perl /opt/vep/src/ensembl-vep/convert_cache.pl --dir ${vep_cache} --species $species --version all
        else
            echo "Downloading VEP cache for species ${species} from ${vep_cache_url}"
            wget -O cache_file.tar.gz ${vep_cache_url}
            tar -xzf cache_file.tar.gz -C ${vep_cache}/
        fi
    fi
    if [ ! -z "${vep_ref_fasta_url}" ]; then
        if [ ! -f "${vep_cache}/${species}/${vep_fasta_path}" ]; then
            echo "Downloading VEP reference fasta for species ${species} from ${vep_ref_fasta_url}"
            wget -O ${vep_cache}/${species}/${vep_fasta_path} ${vep_ref_fasta_url}
            gunzip -c ${vep_cache}/${species}/${vep_fasta_path} | bgzip --index --index-name ${vep_cache}/${species}/${vep_fasta_path}.gzi > ${vep_cache}/${species}/${vep_fasta_path}_2
            mv ${vep_cache}/${species}/${vep_fasta_path}_2 ${vep_cache}/${species}/${vep_fasta_path}
        fi
    fi
    
    # if input_ref_fasta was given, use it as ref_fasta as-is; otherwise (and if default_ref_fasta_url
    # is not empty) download the default reference fasta and use it as ref_fasta
    if [ ! -z "${input_ref_fasta}" ] ; then
        echo "Using the supplied reference fasta ${input_ref_fasta}"
        if [ "${input_ref_fasta}" != "ref_fasta.fa" ] ; then
            ln -s ${input_ref_fasta} ref_fasta.fa
        fi
        if [ ! -z "${input_ref_fasta_index}" ] ; then
            echo "Using the index found next to it: ${input_ref_fasta_index}"
            if [ "${input_ref_fasta_index}" != "ref_fasta.fa.fai" ] ; then
                ln -s ${input_ref_fasta_index} ref_fasta.fa.fai
            fi
        else
            echo "No ${input_ref_fasta}.fai next to the supplied reference fasta, indexing it"
            # htslib's faidx, via the same perl bindings VEP itself uses to read --fasta
            # (the ensembl-vep image ships no samtools); builds ref_fasta.fa.fai in place
            perl -e 'use Bio::DB::HTS::Faidx; Bio::DB::HTS::Faidx->new("ref_fasta.fa");'
        fi
    elif [ ! -f "${vep_cache}/ref_fasta.fa" ] ; then
        if [ ! -z "${default_ref_fasta_url}" ] ; then
            echo "No input reference fasta given, downloading default reference fasta from ${default_ref_fasta_url}"
            wget -O ref_fasta.fa ${default_ref_fasta_url}
            wget -O ref_fasta.fa.fai ${default_ref_fasta_url}.fai
            cp ref_fasta.fa ${vep_cache}/ref_fasta.fa
            cp ref_fasta.fa.fai ${vep_cache}/ref_fasta.fa.fai
        else
            echo "No input reference fasta given and default reference fasta URL is also empty, cannot proceed"
            exit 1
        fi
    else
        echo "No input reference fasta given, using existing reference fasta from VEP cache"
        cp ${vep_cache}/ref_fasta.fa ref_fasta.fa
        cp ${vep_cache}/ref_fasta.fa.fai ref_fasta.fa.fai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
