/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nf-prepare-vcf_pipeline'
include { COMPUTE_GENE_RANGES       } from '../modules/local/python/compute_gene_ranges'
include { BCFTOOLS_ANNOTATE         } from '../modules/local/bcftools/annotate'
include { BCFTOOLS_NORM             } from '../modules/local/bcftools/norm'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1  } from '../modules/local/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2  } from '../modules/local/bcftools/index'
include { VEP_ANNOTATE             } from '../modules/local/vep/annotate'
include { VEP_UPDATECACHE          } from '../modules/local/vep/updatecache'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NF_PREPARE_VCF {

    take:
    ch_input_vcf // channel: VCF file read in from --input_vcf
    
    main:

    ch_versions = Channel.empty()
    vep_cachedir = "${projectDir}/../vep_cachedir"
    ch_vep_cachedir = Channel.fromPath(vep_cachedir, checkIfExists: true)
    ch_rename_chr = Channel.fromPath("${projectDir}/assets/rename_chr.txt", checkIfExists: true)
    ch_meta = ch_input_vcf.map { t -> t[0] }
    

    VEP_UPDATECACHE (
        ch_meta.first().combine(ch_vep_cachedir),
        Channel.value(params.vep_updatecache_species),
        Channel.value(params.vep_updatecache_options),
        Channel.value(params.vep_cache_url),
        Channel.value(tuple(params.ref_fasta_url, params.vep_fasta_path))
    )
    ch_vep_cachesubdir = VEP_UPDATECACHE.out.cachesubdir.first()
    ch_versions = ch_versions.mix(VEP_UPDATECACHE.out.versions.first())


    BCFTOOLS_INDEX_1 (
        ch_input_vcf
    )
    ch_input_vcf_tbi = BCFTOOLS_INDEX_1.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions.first())


    BCFTOOLS_ANNOTATE (
        ch_input_vcf
            .join(ch_input_vcf_tbi, by: 0)  // Join by the first element (meta)
            .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file, [], [], [], []) }
            .combine(ch_rename_chr),
        Channel.value("rename_chr"),
    )
    ch_annotated_vcf = BCFTOOLS_ANNOTATE.out.vcf
    ch_annotated_vcf_tbi = BCFTOOLS_ANNOTATE.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    BCFTOOLS_NORM (
        ch_annotated_vcf
            .join(ch_annotated_vcf_tbi, by: 0)
            .combine(ch_vep_cachesubdir.map { t -> "${t}/${params.vep_fasta_path}" })
            .map { meta, vcf_file, tbi_file, fasta_file -> tuple(meta, vcf_file, tbi_file, fasta_file, []) },
        Channel.value("norm"),
    )
    ch_normalized_vcf = BCFTOOLS_NORM.out.vcf
    ch_normalized_vcf_tbi = BCFTOOLS_NORM.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())
    ch_tracking = BCFTOOLS_NORM.out.tracking_out.first()


    VEP_ANNOTATE (
        ch_normalized_vcf
            .join(ch_normalized_vcf_tbi, by: 0)
            .combine(ch_vep_cachesubdir),
        Channel.value(params.vep_annotate_species),
        Channel.value(params.vep_fasta_path),
        Channel.value(params.vep_annotate_options)
    )
    ch_vep_vcf  = VEP_ANNOTATE.out.vcf
    ch_versions = ch_versions.mix(VEP_ANNOTATE.out.versions.first())


    BCFTOOLS_INDEX_2 (
        ch_vep_vcf
    )
    ch_vep_vcf_tbi = BCFTOOLS_INDEX_2.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())


    compute_ranges_py_script_ch = Channel.fromPath(params.compute_ranges_py_script_path, checkIfExists: true)
    COMPUTE_GENE_RANGES (
        ch_vep_vcf
            .combine(compute_ranges_py_script_ch),
        Channel.value("compute_ranges")
    )
    ch_all_gene_ranges = COMPUTE_GENE_RANGES.out.gene_ranges
    ch_versions = ch_versions.mix(COMPUTE_GENE_RANGES.out.versions.first())

    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nf-prepare-vcf_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
