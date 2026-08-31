# nf-prepare-vcf: Output

Everything the pipeline writes, and how to read it. All paths are relative to
`--outdir`.

Most directories appear only when `--publish_intermediate true` is set. The tables below
say which is which.

## Written by default

### `bcftools_reheader/` and `bcftools_index` - the prepared VCF and its index

This is the final VCF after all the preparation steps.

### `compute_gene_ranges/` - the gene ranges calculated from the VCF

The gene ranges are used by [nf-eval-gene-assoc](https://github.com/biodatageeks/nf-eval-gene-assoc) and are not needed
for [nf-rare-var-assoc](https://github.com/biodatageeks/nf-rare-var-assoc) operation.

### `plink2_ld_report/` - Linkage Disequilibrium report for the VCF

The LD report is an optional input of [nf-eval-gene-assoc](https://github.com/biodatageeks/nf-eval-gene-assoc) and is
not needed for [nf-rare-var-assoc](https://github.com/biodatageeks/nf-rare-var-assoc) operation.

Note: the `compute_ld_report` flag is by default `false`, turning it on will trigger the LD report computation, which
may take a lot of disc space.

### `pipeline_info/`

Nextflow's own records: `execution_report_*.html` (per-task runtime and memory),
`execution_timeline_*.html`, `execution_trace_*.txt` and `pipeline_dag_*.html`. Use the
trace file to find which steps are slow or run out of memory.

## Written only with `--publish_intermediate true`

These are intermediate results, not needed during normal usage of nf-prepare-vcf. For a list of modules producing
those results please see `conf/modules.config` and look for `enabled: params.publish_intermediate`.