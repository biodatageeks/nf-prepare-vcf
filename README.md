# nf-prepare-vcf

[![GitHub Actions CI Status](https://github.com/biodatageeks/nf-prepare-vcf/actions/workflows/nf-test.yml/badge.svg)](https://github.com/biodatageeks/nf-prepare-vcf/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/biodatageeks/nf-prepare-vcf/actions/workflows/linting.yml/badge.svg)](https://github.com/biodatageeks/nf-prepare-vcf/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**nf-prepare-vcf** is a Nextflow pipeline that is a sub-pipeline of the **nf-rare-var-assoc** rare-variant association
studies pipeline. nf-prepare-vcf ingests a multi-sample VCF file and executes data preprocessing and preparation steps:

- **Sorting, splitting multi-allelic sites and removing exact duplicates.** A site with
  several alternative alleles becomes one record per allele, so that each can be tested
  and annotated separately.
- **Chromosome renaming and variant identifiers.** Every variant is given the identifier
  `CHROM_POS_REF_ALT`. This happens *after* the multi-allelic split, so that identifiers
  stay unique. All later steps match variants to genes by this identifier, so it has to
  be assigned consistently.
- **Left-alignment is deliberately switched off** (`--do-not-normalize`). Left-aligning
  indels would shift their positions and therefore change their identifiers, breaking
  the correspondence with the annotation files. This is a deliberate departure from the
  usual convention, made so that variant identifiers stay stable through the pipeline.
- **VEP annotation.** Predicted consequences are written into the VCF's `CSQ` field, and
  are what the gene grouping step later reads.
- **Dosage computation from genotype likelihoods.** A hard genotype call throws away how
  confident the caller was. Instead, a dosage - the expected number of alternative
  alleles, between 0 and 2 - is computed from the `PL` genotype likelihoods and stored
  in a `DS` field. The computation also corrects calls where every likelihood is zero,
  which some callers emit for homozygous-reference sites, using the genotype quality
  instead. Only genotypes above a minimum quality contribute.
- **Gene ranges computation.**

## Usage

> [!NOTE]
> See [Usage](docs/usage.md) and [Output](output.md) for more details.

```bash
nextflow run biodatageeks/nf-prepare-vcf \
   -profile <docker/singularity/...> \
   --input_vcf input.vcf.gz \
   --project_name myproject \
   --outdir <OUTDIR>
```

## Citations

This pipeline uses code and infrastructure developed and maintained by the
[nf-core](https://nf-co.re) community, reused here under the
[MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

## License

Copyright (c) 2026 Piotr Suszyński.

This pipeline is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. See the
[LICENSE](LICENSE) file for the full text.

It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.