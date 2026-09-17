# `tests/fixtures/`

Integration-test fixtures for `nf-prepare-vcf`.

## `unprepared_rand_500.vcf.gz`

Raw, *unprepared* 1000 Genomes high-coverage VCF used by `tests/prepare_invariants.nf.test`
to exercise the full preparation pipeline on messy input. It is also the `-profile test` input
(`conf/test.config`): it replaced the old sim_chr22 dataset, which was already clean --
biallelic, single chromosome, IDs preset -- so it could not prove the NORM/ANNOTATE
transformations actually fire, and which was not `chr`-prefixed as BCFTOOLS_NORM now requires.

- 3 chromosomes incl. X: `chr12`, `chr22`, `chrX` (with the `chr` prefix)
- 3202 samples (1000G high-coverage release, IDs `HG00096`..`NA21144`)
- ~500 variants, including 35 multiallelic sites and duplicate variant IDs
- carries per-genotype `PL` (used to validate `CALC_DOSAGE_POLARSBIO`'s `DS` against spec)

## `deepvariant_real_2samples.vcf.gz`

The only real DeepVariant output available for this work: 11 records on `chr22`, two samples,
used by the dosage module suite to check caller detection and the flat prior on a header no one
crafted for the test. It carries `VAF`, `MED_DP` and `MIN_DP` — the last of which is why `MIN_DP`
is not part of the GATK fingerprint in `calc_dosage.py` — and `##DeepVariant_version=1.8.0`, but
no `##source` line, so detection has to fall through to the `FORMAT` fingerprint.

Built from nf-core test-dataset output:

```bash
DV=/data/doktorat/biodatageeks/nf_out_1__2bams_nf_core_test_datasets/deepvariant
bcftools merge -Ou $DV/test1.vcf.gz $DV/test2.vcf.gz \
  | bcftools +fill-tags -Oz -o deepvariant_real_2samples.vcf.gz -- --tags AC,AN
bcftools index -t deepvariant_real_2samples.vcf.gz
```

The merge is not cosmetic: polars-bio reads a single-sample VCF into flat `GT`/`PL` columns
instead of the `genotypes` struct the dosage query is written against, so `calc_dosage.py`
cannot read one at all. `+fill-tags` supplies the `INFO/AC`/`AN` the prior needs, exactly as
`BCFTOOLS_NORM` does in the pipeline.

## Notes on `unprepared_rand_500.vcf.gz`

This file is copied verbatim from the sibling `nf-rare-var-assoc` repo
(`assets/three_chr_unprepared/unprepared_rand_500.vcf.gz`), where it is the canonical
unprepared test fixture. It was brought here during PB5 of that repo's test-cleanup plan so
`nf-prepare-vcf` has self-contained integration coverage of its preparation invariants instead
of relying on the downstream repo's wrapper test.
