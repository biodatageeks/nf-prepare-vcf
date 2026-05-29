# `tests/fixtures/`

Integration-test fixtures for `nf-prepare-vcf`.

## `unprepared_rand_500.vcf.gz`

Raw, *unprepared* 1000 Genomes high-coverage VCF used by `tests/prepare_invariants.nf.test`
to exercise the full preparation pipeline on messy input (the `-profile test` sim_chr22
dataset is already clean -- biallelic, single chromosome, IDs preset -- so it cannot prove the
NORM/ANNOTATE transformations actually fire).

- 3 chromosomes incl. X: `chr12`, `chr22`, `chrX` (with the `chr` prefix)
- 3202 samples (1000G high-coverage release, IDs `HG00096`..`NA21144`)
- ~500 variants, including 35 multiallelic sites and duplicate variant IDs
- carries per-genotype `PL` (used to validate `CALC_DOSAGE_POLARSBIO`'s `DS` against spec)

This file is copied verbatim from the sibling `nf-rare-var-assoc` repo
(`assets/three_chr_unprepared/unprepared_rand_500.vcf.gz`), where it is the canonical
unprepared test fixture. It was brought here during PB5 of that repo's test-cleanup plan so
`nf-prepare-vcf` has self-contained integration coverage of its preparation invariants instead
of relying on the downstream repo's wrapper test.
