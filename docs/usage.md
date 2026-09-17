# nf-prepare-vcf: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

You can run the pipeline using:

```bash
nextflow run biodatageeks/nf-prepare-vcf \
   -profile <docker/singularity/...> \
   --input_vcf input.vcf.gz \
   --project_name myproject \
   --outdir <OUTDIR>
```

The input VCF file, optionally compressed, should be a multi-sample VCF - it can be obtained from gvcf files for example by
running the [GLNexus](https://github.com/dnanexus-rnd/GLnexus) joint variant genotyping software. The VCF file should also
contain an allele frequency field in the INFO column.

For genotype dosages to be computed, the FORMAT field must contain genotype likelihoods (`PL`) alongside `GT`, `GQ` and `DP`.
Without them the pipeline still runs, but on hard genotype calls only (in such case `GT`, `GQ` and `DP` are mandatory).

### Genotype dosage options

A dosage is the expected number of alternative alleles, computed from `PL` by Bayes' theorem. What
the pipeline assumes `PL` means decides which prior it applies, and that is the single largest
influence on the resulting dosages.

- `--variant_caller` (default: empty) — the caller that produced the input, used to decide how to
  read `PL`. Empty means detect it from the input header: a caller named in `##source` or
  `##GATKCommandLine`, otherwise `FORMAT` fields that only one caller emits. A file that identifies
  no caller is assumed to be GATK-like, and the pipeline says so loudly in the log.
- `--pl_type` (default: `auto`) — overrides the above. `likelihood` (GATK-like) means `PL` holds
  pure genotype likelihoods, as the VCF specification defines them, so a Hardy-Weinberg prior at
  the variant's allele frequency is applied. `posterior` (DeepVariant-like) means `PL` already
  carries the caller's own prior, so it is only normalised.
- `--ds_prior_min_af` (default: 0) — floor on the allele frequency used in the prior. It softens
  how strongly a rare variant pulls a weak carrier's dosage towards 0. The estimate is already
  floored at `1/2 / (AN + 1)`, so 0 does not mean a zero frequency.
- `--ds_max_pl_min` (default: 30) — genotypes whose smallest `PL` value reaches this get no dosage.
  Splitting a multiallelic record leaves a sample that carries some *other* alternative allele
  with a `PL` that has no zero in it, because the entry that was 0 belonged to a genotype the
  split record cannot express. Such a sample supports neither this alternative allele nor the
  reference, so it is left out rather than dosed; below the threshold its dosage is weighted
  between the likelihoods and its `GT`. Raise it to dose more of these genotypes, lower it to
  drop them sooner.
- `--force_dosage_recalc` (default: `false`) — by default a `FORMAT/DS` already declared in the
  input header is trusted and passed through untouched. This recomputes every dosage instead.
- `--calc_ds_min_gq` (default: 3) — genotypes below this `GQ` get no dosage.

Haploid genotypes — male chrX outside the pseudo-autosomal regions, chrY, the mitochondrial
genome — are dosed on the haploid scale, so their `DS` lies in `[0, 1]` rather than `[0, 2]`. That
is the VCF reading of `DS` as the expected number of alternative alleles, and it is what keeps the
dosage on the same scale as the hard call beside it: PLINK 2 reads the ploidy off `GT`, not off
the sample's declared sex, and doubles a haploid dosage on import, exactly as it counts a
hemizygous hard call as two copies.

The dosage step writes a tracking JSON beside its output, as `BCFTOOLS_NORM` does, saying how many
genotypes were computed, how many were reconstructed from `GQ`, and how many were left without a
dosage, along with any warning the run raised. When the input's own dosages are trusted the query
never runs, so the counts are all `null` and `passthrough` is `true`. It is published with the
other intermediates, under `--publish_intermediate`.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
