# Test fixtures for `bcftools/norm`

All fixtures use chromosome **MT192765.1** (sarscov2) so they can be used with the nf-core
sarscov2 genome FASTA (`genomics/sarscov2/genome/genome.fasta.gz`) as `--fasta-ref`.
REF alleles were extracted from that FASTA at the respective positions.

| File | Records | Description |
|---|---|---|
| `multiallelic_sarscov2.vcf.gz` | 8 | 5 biallelic + 3 multiallelic sites (pos 200: 2 ALTs, 400: 2 ALTs, 600: 3 ALTs). After `bcftools norm -m -any` → 12 records. |
| `split_sarscov2.vcf.gz` | 12 | All-biallelic, pre-split form of the above. After `bcftools norm -m +any` → 8 records. |
| `multiallelic_samples_sarscov2.vcf.gz` | 8 | The same 8 records with three diploid samples added, so the `bcftools +fill-tags` the module pipes into has genotypes to count `AC`/`AN` from. The other two carry no samples at all, where `+fill-tags` writes no tags. |

Genotypes of `multiallelic_samples_sarscov2.vcf.gz`, and the `AC`/`AN` they imply after the
split, which is what the test asserts:

| POS | REF | ALT | S1 | S2 | S3 | after `-m -any` |
|---|---|---|---|---|---|---|
| 100 | A | T | 0/0 | 0/1 | 1/1 | `AC=3;AN=6` |
| 200 | C | T,G | 0/1 | 0/2 | 1/2 | T: `AC=2`, G: `AC=2` |
| 300 | C | A | 0/0 | 0/0 | 0/1 | `AC=1` |
| 400 | G | C,A | 0/1 | 0/0 | 0/2 | C: `AC=1`, A: `AC=1` |
| 500 | A | G | 1/1 | 1/1 | 1/1 | `AC=6` |
| 600 | A | T,G,C | 0/1 | 0/2 | 0/3 | T, G, C: `AC=1` each |
| 700 | G | A | ./. | 0/1 | 0/1 | `AC=2;AN=4` — a missing genotype lowers `AN` |
| 800 | C | T | 0/0 | 0/0 | 0/0 | `AC=0` |
