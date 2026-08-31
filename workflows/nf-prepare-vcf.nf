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
include { BCFTOOLS_SORT             } from '../modules/local/bcftools/sort'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1  } from '../modules/local/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2  } from '../modules/local/bcftools/index'
include { BCFTOOLS_REHEADER         } from '../modules/local/bcftools/reheader'
include { BCFTOOLS_VCF2PSAM         } from '../modules/local/bcftools/vcf2psam'
include { PLINK2_MAKEPGEN           } from '../modules/local/plink2/makepgen'
include { PLINK2_LD_REPORT          } from '../modules/local/plink2/ld_report'
include { VEP_ANNOTATE              } from '../modules/local/vep/annotate'
include { VEP_UPDATECACHE           } from '../modules/local/vep/updatecache'
include { CALC_DOSAGE_POLARSBIO     } from '../modules/local/python/calc_dosage_polarsbio'



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
    python_calc_dosage_polarsbio_script_ch = Channel.fromPath(params.calc_dosage_polarsbio_script_path, checkIfExists: true)
    

    VEP_UPDATECACHE (
        ch_meta.first().combine(ch_vep_cachedir),
        Channel.value(params.vep_updatecache_species),
        Channel.value(params.vep_updatecache_options),
        Channel.value(params.vep_cache_url),
        Channel.value(tuple(params.ref_fasta_url, params.vep_fasta_path))
    )
    ch_vep_cachesubdir = VEP_UPDATECACHE.out.cachesubdir.first()
    ch_versions = ch_versions.mix(VEP_UPDATECACHE.out.versions.first())


    BCFTOOLS_SORT (
        ch_input_vcf,
        Channel.value("sort"),
    )
    ch_sorted_vcf = BCFTOOLS_SORT.out.vcf
    ch_sorted_vcf_tbi = BCFTOOLS_SORT.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions.first())


    // NORM must be before ANNOTATE (where we assign variant ids) because we must first split multiallelic sites before assigning variant ids.
    // Otherwise we'll end up with duplicated ids with a comma within them and this will cause subsequent plink write-snplist steps to fail
    BCFTOOLS_NORM (
        ch_sorted_vcf
            .join(ch_sorted_vcf_tbi, by: 0)
            .combine(ch_vep_cachesubdir.map { t -> "${t}/${params.vep_fasta_path}" })
            .map { meta, vcf_file, tbi_file, fasta_file -> tuple(meta, vcf_file, tbi_file, fasta_file, []) },
        Channel.value("norm"),
    )
    ch_normalized_vcf = BCFTOOLS_NORM.out.vcf
    ch_normalized_vcf_tbi = BCFTOOLS_NORM.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())
    ch_tracking = BCFTOOLS_NORM.out.tracking_out.first()
    
    BCFTOOLS_ANNOTATE (
        ch_normalized_vcf
            .join(ch_normalized_vcf_tbi, by: 0)  // Join by the first element (meta)
            .map { meta, vcf_file, tbi_file -> tuple(meta, vcf_file, tbi_file, [], [], [], []) }
            .combine(ch_rename_chr),
        Channel.value("rename_chr"),
    )
    ch_annotated_vcf = BCFTOOLS_ANNOTATE.out.vcf
    ch_annotated_vcf_tbi = BCFTOOLS_ANNOTATE.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())


    VEP_ANNOTATE (
        ch_annotated_vcf
            .join(ch_annotated_vcf_tbi, by: 0)
            .combine(ch_vep_cachesubdir),
        Channel.value(params.vep_annotate_species),
        Channel.value(params.vep_fasta_path),
        Channel.value(params.vep_annotate_options)
    )
    ch_vep_vcf  = VEP_ANNOTATE.out.vcf
    ch_versions = ch_versions.mix(VEP_ANNOTATE.out.versions.first())


    CALC_DOSAGE_POLARSBIO (
        ch_vep_vcf.combine(python_calc_dosage_polarsbio_script_ch),
        Channel.value(params.calc_ds_min_gq),
        Channel.value("calc_dosage")
    )
    ch_with_dosage_vcf = CALC_DOSAGE_POLARSBIO.out.vcf
    ch_versions = ch_versions.mix(CALC_DOSAGE_POLARSBIO.out.versions.first())


    BCFTOOLS_REHEADER (
        ch_vep_vcf.join(ch_with_dosage_vcf, by: 0),
        Channel.value("##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage (expected number of non-reference alleles)\">"),
        Channel.value("reheader")
    )
    ch_with_correct_header_vcf = BCFTOOLS_REHEADER.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions.first())


    BCFTOOLS_INDEX_2 (
        ch_with_correct_header_vcf
    )
    ch_with_correct_header_vcf_tbi = BCFTOOLS_INDEX_2.out.tbi
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())


    compute_ranges_py_script_ch = Channel.fromPath(params.compute_ranges_py_script_path, checkIfExists: true)
    COMPUTE_GENE_RANGES (
        ch_with_correct_header_vcf
            .combine(compute_ranges_py_script_ch),
        Channel.value("compute_ranges")
    )
    ch_all_gene_ranges = COMPUTE_GENE_RANGES.out.gene_ranges
    ch_versions = ch_versions.mix(COMPUTE_GENE_RANGES.out.versions.first())


    BCFTOOLS_VCF2PSAM (
        ch_with_correct_header_vcf.join(ch_with_correct_header_vcf_tbi, by: 0)
    )
    ch_unk_sex_psam = BCFTOOLS_VCF2PSAM.out.psam
    ch_versions = ch_versions.mix(BCFTOOLS_VCF2PSAM.out.versions.first())


    PLINK2_MAKEPGEN (
        ch_with_correct_header_vcf
            .join(ch_with_correct_header_vcf_tbi, by: 0)
            .join(ch_unk_sex_psam, by: 0)
            .map { meta, vcf_file, tbi_file, unk_sex_psam_file -> tuple(meta, [], [], unk_sex_psam_file, vcf_file, tbi_file, [], [], [], []) },
        Channel.value(''),
        Channel.value(''),
        Channel.value(params.plink2_makepgen_vcf_input_options),
        Channel.value('makepgen'),
        Channel.value(params.plink2_makepgen_options)
    )
    ch_pgen_pvar_psam  = PLINK2_MAKEPGEN.out.out_pgen_pvar_psam
    ch_versions = ch_versions.mix(PLINK2_MAKEPGEN.out.versions.first())
    ch_tracking = ch_tracking.mix(PLINK2_MAKEPGEN.out.tracking_out.first())


    if (params.compute_ld_report == true) {
        PLINK2_LD_REPORT (
            ch_pgen_pvar_psam
                .map { meta, pgen, pvar, psam -> tuple(meta, pgen, pvar, psam, []) },
            Channel.value('ld_report'),
            Channel.value(params.plink2_ld_report_options)
        )
        ch_ld_report  = PLINK2_LD_REPORT.out.ld_report
        ch_versions = ch_versions.mix(PLINK2_LD_REPORT.out.versions.first())
    }


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
