"""Compute FORMAT/DS, the expected number of alternative alleles, for every genotype.

The maths is written below in the usual symbols; the code uses these names for them:
m -> ploidy, g -> the genotype as a count of alt alleles, L_g -> likelihood_ref/_het/_alt,
pi_g -> prior_ref/_het/_alt, pi_g*L_g -> post_ref/_het/_alt, p -> allele_freq,
w -> call_err_prob, DS -> ds.

PL is Phred-scaled and normalised, so L_g = 10^(-PL_g/10) is the likelihood of genotype g up
to a constant shared by all genotypes of that sample. Bayes' theorem with a prior pi_g gives
the posterior P(g|D) = pi_g L_g / sum_h pi_h L_h, and the dosage is its mean:

    DS = sum_g g pi_g L_g / sum_g pi_g L_g     (diploid g in {0,1,2}, haploid g in {0,1})

The prior is what --pl-type selects, and it is the single largest influence on the result:

  * likelihood (GATK-like, the default): PL holds pure genotype likelihoods, which is what
    the VCF specification defines it to be, so the population prior still has to be supplied.
    It is Hardy-Weinberg at the allele frequency p estimated per record from INFO/AC and
    INFO/AN: diploid pi = ((1-p)^2, 2p(1-p), p^2), haploid pi = (1-p, p).
  * posterior (DeepVariant-like): PL is derived from the caller's own posterior
    probabilities, which already carry a prior. Applying a second one would count the
    population information twice, so the prior is flat and the formula above reduces to
    normalising the likelihoods.

Under --pl-type auto, the default, the caller decides: --variant-caller if it was given, else a
caller named in the header's provenance lines, else FORMAT fields that only one caller emits.
A file that identifies no caller is assumed to be GATK-like, with a warning.

Reconstruction from GQ for GT=0/0. A genotype called homozygous reference whose PL is all zeros
or absent says nothing about *which* genotype the truth is, so GQ is all that is left. Under
the VCF definition of GQ, w = 10^(-GQ/10) is the probability that the call is wrong, i.e.
w = P(0/1|D) + P(1/1|D). Nothing informs how w divides between those two, so it is divided in
proportion to their priors - which genotype the truth is, is a population question, so this
split is Hardy-Weinberg whatever --pl-type says:

    P(0/1|D) = w * pi_1/(pi_1+pi_2) = w * 2p(1-p)/(2p(1-p) + p^2) = w * 2(1-p)/(2-p)
    P(1/1|D) = w * pi_2/(pi_1+pi_2) = w *     p^2/(2p(1-p) + p^2) = w *      p/(2-p)
    DS       = P(0/1|D) + 2 P(1/1|D) = w * (2(1-p) + 2p)/(2-p)    = 2w/(2-p)

Both priors share a factor p, which cancels; that is why the frequency survives only in the
single term 2-p, and why DS stays close to w for any p. For a haploid call the whole of w
goes to the one other genotype, so DS = w. w is capped at m/(m+1), because a caller reports
the most probable genotype and that genotype therefore keeps at least 1/(m+1).

PL without a zero, after a multiallelic split. BCFTOOLS_NORM splits a multiallelic record into
one record per ALT, keeping for each the PL entries of the m+1 genotypes over REF and that ALT,
and rewriting the other ALTs in GT as REF. A sample whose own genotype involves a different ALT
therefore keeps none of its PL zero - the entry that was 0 belonged to a genotype this record
cannot express - which breaks the VCF convention that the most likely genotype has PL 0. The
likelihoods that remain are real but they only account for part of the sample's probability
mass: L_min = 10^(-min(PL)/10) of it, the rest sitting on genotypes the record cannot express.
So the two parts are weighted (pl_min, likelihood_weight, gt_alt_count, ds_value in the query):

    DS = L_min * DS_likelihood + (1 - L_min) * (copies of this ALT in GT)

The second term is what GT contributes, since only it says what those other genotypes mean for
this ALT: bcftools resolved that when it split the record. A well-formed PL has min(PL) = 0,
hence L_min = 1, and the formula is the plain dosage above. Once min(PL) reaches
--ds-max-pl-min, DS is left missing instead: such a sample supports neither this ALT nor the
reference - its GT reads 0/0 only because the split had nowhere else to put it - and any number
would tell a downstream test that it is a confident reference carrier.

For example, take a site with four alleles, REF=A and ALT=C,G,T, and a sample whose genotype is
C/G. Its PL has ten entries, one per unordered genotype, and the 0 sits on C/G. The split gives
three biallelic records, each keeping only the three entries over REF and its own ALT:

    record ALT=C   GT 1/0   PL 300,40,12     min 12    truth: 1 copy of C
    record ALT=G   GT 0/1   PL 300,44,14     min 14    truth: 1 copy of G
    record ALT=T   GT 0/0   PL 300,260,247   min 247   truth: 0 copies of T

Not one of them keeps a zero, because C/G is gone from all three. In the first record the
likelihoods on their own read as C/C - 12 is the smallest of them, since a sample with no A at
all is punished hardest by the genotypes containing one - and the dosage comes out near 2; the
weighting pulls it back to about 1, which is what GT says. The last record is the same effect
at its extreme: T/T is the least bad of three bad options, so the likelihoods alone would dose
this sample near 2 for an allele it does not carry. GT alone would say 0, but the sample has no
A either, so it is not a reference carrier to be compared against T carriers. DS is left
missing.
"""

import argparse
import gzip
import json
import shutil
from pathlib import Path
import polars as pl
import polars_bio as pb


# What each known caller means for PL. DeepVariant comes first because it is matched against the
# provenance lines as a substring: a DeepVariant call set later handled by a GATK tool carries both
# names, and it is the caller that decides what PL holds, not whatever ran afterwards.
PL_TYPE_BY_CALLER = {
    "deepvariant": "posterior",
    "gatk": "likelihood",
    "haplotypecaller": "likelihood",
}

# FORMAT fields only one of the two emits. This is the fallback fingerprint for the common case of
# a header that names no caller at all - the pipeline's own reference file is such a file.
CALLER_FORMAT_FINGERPRINTS = [
    ("DeepVariant", "posterior", {"VAF", "MED_DP"}),
    ("GATK", "likelihood", {"RGQ", "PGT", "PID"}),
]

parser = argparse.ArgumentParser()
parser.add_argument("--input-vcf-path", help="Input VCF file")
parser.add_argument("--output-vcf-path", help="Output VCF file")
parser.add_argument("--calc-ds-min-gq", type=int, default=1, help="Minimum genotype quality for DS calculation")
parser.add_argument(
    "--pl-type",
    choices=["likelihood", "posterior", "auto"],
    default="auto",
    help=(
        "How to read PL. 'likelihood' (GATK-like, and what the VCF specification defines PL to "
        "be) means PL holds pure genotype likelihoods, so a Hardy-Weinberg prior is applied. "
        "'posterior' (DeepVariant-like) means PL already carries the caller's own prior, so it is "
        "only normalised. 'auto' decides from the caller named in the input header, falling back "
        "to 'likelihood'."
    ),
)
parser.add_argument(
    "--ds-prior-min-af",
    type=float,
    default=0.0,
    help=(
        "Lower bound on the allele frequency used in the Hardy-Weinberg prior. Note that 0 does "
        "not allow p = 0: the estimate (AC + 1/2) / (AN + 1) already floors itself at "
        "1/2 / (AN + 1)."
    ),
)
parser.add_argument(
    "--ds-max-pl-min",
    type=int,
    default=30,
    help=(
        "Leave DS missing once the smallest PL value of a genotype reaches this, i.e. once the "
        "sample's own genotype is at least 10^(this/10) times more likely than anything the "
        "record can express. Only a multiallelic split produces such a genotype; see the module "
        "docstring."
    ),
)
parser.add_argument(
    "--variant-caller",
    default="",
    help=(
        "Name of the variant caller that produced the input, used to resolve --pl-type auto "
        f"without inspecting the header. One of: {', '.join(PL_TYPE_BY_CALLER)}. Empty (the "
        "default) means detect it from the header."
    ),
)
parser.add_argument(
    "--force-dosage-recalc",
    action="store_true",
    help=(
        "Compute a dosage for every genotype even when the input file already declares FORMAT/DS, "
        "overwriting the input values."
    ),
)
parser.add_argument(
    "--tracking-json-path",
    default="",
    help="Where to write the tracking JSON: how many genotypes ended up in each case, and any "
         "warning raised. Empty (the default) writes none.",
)
parser.add_argument(
    "--process-name",
    default="CALC_DOSAGE_POLARSBIO",
    help="Name recorded in the tracking JSON, so it reads like the other steps' tracking files.",
)
parser.add_argument("--workflow-name", default="", help="Workflow recorded in the tracking JSON.")

args = parser.parse_args()

warnings: list[str] = []

# What became of each genotype. "passed_through" is only ever non-zero when the input's own
# dosages are kept, and then nothing is counted at all, so it is there for the shape, not for a
# value; "reconstructed" is the GQ branch, for a hom-ref call whose PL says nothing.
GENOTYPE_COUNTS = ["total", "passed_through", "computed", "reconstructed", "missing"]


def warn(message: str) -> None:
    """Warnings go to the tracking file as well as the log, so a test can read them from there."""
    warnings.append(message)
    print(f"WARNING: {message}")


def write_tracking(passthrough: bool) -> None:
    """The counts of what happened to each genotype, in the shape BCFTOOLS_NORM's tracking uses.

    Reporting is not worth failing a run for: by the time this is called the output VCF is
    written, so anything that goes wrong here is logged and swallowed.
    """
    if not args.tracking_json_path:
        return
    try:
        # In passthrough mode the input's own dosages were kept, the query never ran, and there
        # is nothing to count: counting would mean a pass over the file for nothing.
        counts = dict.fromkeys(GENOTYPE_COUNTS, None)
        if not passthrough:
            counted = pb.sql(COUNTS_SQL).collect().to_dicts()[0]
            counts = {name: int(counted[name]) for name in GENOTYPE_COUNTS}
            print(f"Genotypes: {counts}")
        Path(args.tracking_json_path).write_text(json.dumps({
            "process_name": args.process_name,
            "workflow_name": args.workflow_name,
            "parameters": (
                f"--calc-ds-min-gq {args.calc_ds_min_gq} --pl-type {args.pl_type} "
                f"--ds-prior-min-af {args.ds_prior_min_af} --ds-max-pl-min {args.ds_max_pl_min}"
                + (f" --variant-caller {args.variant_caller}" if args.variant_caller else "")
                + (" --force-dosage-recalc" if args.force_dosage_recalc else "")
            ),
            "passthrough": passthrough,
            "genotypes": counts,
            "warnings": warnings,
        }, indent=2) + "\n")
        print(f"Wrote: {args.tracking_json_path}")
    except Exception as failure:
        print(
            f"WARNING: could not write the tracking file {args.tracking_json_path}: {failure!r}. "
            "The dosages themselves are unaffected."
        )

if not 0.0 <= args.ds_prior_min_af < 1.0:
    parser.error(f"--ds-prior-min-af must lie in [0, 1), got {args.ds_prior_min_af}")

if args.ds_max_pl_min < 1:
    # 0 would drop every genotype: a well-formed PL has a smallest value of exactly 0.
    parser.error(f"--ds-max-pl-min must be at least 1, got {args.ds_max_pl_min}")


# Read before register_vcf: a FORMAT/DS in the input header decides whether the query runs at all.
source_lf = pb.scan_vcf(args.input_vcf_path)
header = pb.get_metadata(source_lf)["header"]
input_format_fields = header["format_fields"]
input_ds = input_format_fields.get("DS")

# Per-file rule, not per-genotype: the input's dosages are trusted whole, or ours are.
if input_ds is not None and not args.force_dosage_recalc:
    ds_number = str(input_ds.get("number", ""))
    if ds_number != "A":
        # Warn, but still trust: the values may be wrong, that is not grounds to overrule them.
        warn(
            f"the input declares FORMAT/DS with Number={ds_number}, not Number=A. "
            "bcftools splits a per-allele FORMAT field correctly only when it is declared "
            "Number=A, so after a multiallelic split these values may not be the dosage of the "
            "individual allele. They are passed through unchanged regardless."
        )
    # Nothing is computed and --calc-ds-min-gq has no effect, so there is nothing to rewrite.
    print(
        "Input header declares FORMAT/DS: trusting the file, copying it through unchanged. "
        "Pass --force-dosage-recalc to compute dosages instead."
    )
    shutil.copyfile(args.input_vcf_path, args.output_vcf_path)
    print(f"Wrote: {args.output_vcf_path}")
    write_tracking(passthrough=True)
    raise SystemExit(0)

if input_ds is not None:
    print("--force-dosage-recalc given: every dosage in the input will be overwritten.")
else:
    print("Input header does not declare FORMAT/DS: computing a dosage for every genotype.")


def read_meta_lines(path: str) -> list[str]:
    """The header's ## lines. get_metadata exposes the structured ones only, and the caller is
    named in free-text lines (##source, ##GATKCommandLine) that it drops."""
    with open(path, "rb") as probe:
        is_gzip = probe.read(2) == b"\x1f\x8b"
    with (gzip.open if is_gzip else open)(path, "rt", errors="replace") as handle:
        return [line.rstrip("\n") for line in handle if line.startswith("##")]


def resolve_pl_type(requested: str, variant_caller: str) -> tuple[str, str]:
    """Decide whether PL holds likelihoods (Hardy-Weinberg prior) or posteriors (flat prior),
    and say what decided it."""
    if requested != "auto":
        return requested, f"--pl-type {requested}"

    if variant_caller:
        named = PL_TYPE_BY_CALLER.get(variant_caller.strip().lower())
        if named is None:
            parser.error(
                f"--variant-caller {variant_caller} is not a caller this script knows. Use one of "
                f"{', '.join(PL_TYPE_BY_CALLER)}, or set --pl-type directly."
            )
        return named, f"--variant-caller {variant_caller}"

    provenance = [
        line.lower() for line in read_meta_lines(args.input_vcf_path)
        if line.lower().startswith(("##source=", "##gatkcommandline"))
    ]
    for caller, caller_pl_type in PL_TYPE_BY_CALLER.items():
        if any(caller in line for line in provenance):
            return caller_pl_type, f"{caller} named in the header"

    declared_format_fields = set(input_format_fields)
    for caller, caller_pl_type, fingerprint in CALLER_FORMAT_FINGERPRINTS:
        found = sorted(fingerprint & declared_format_fields)
        if found:
            return caller_pl_type, f"the {caller}-specific FORMAT field(s) {', '.join(found)}"

    warn(
        "--pl-type auto could not identify the variant caller - the header names none "
        "and declares no FORMAT field characteristic of one - so PL is assumed to hold "
        "likelihoods (GATK-like, the VCF-specification reading) and a Hardy-Weinberg prior is "
        "used. Pass --pl-type posterior, or --variant-caller deepvariant, for DeepVariant-like "
        "input whose PL already carries a prior: the prior is the single largest influence on "
        "the dosages."
    )
    return "likelihood", "the default assumption, no caller identified"


pl_type, pl_type_reason = resolve_pl_type(args.pl_type, args.variant_caller)
print(
    f"Using the {'Hardy-Weinberg' if pl_type == 'likelihood' else 'flat'} prior "
    f"(--pl-type {pl_type}, from {pl_type_reason})."
)


pb.register_vcf(
    args.input_vcf_path,
    name="vcf_table",
    # info_fields unset: registers every INFO field, which is what keeps them in the output.
    format_fields=["GT", "DP", "GQ", "PL"],
)

schema_df = pb.sql("DESCRIBE vcf_table").collect()
name_col = next((c for c in ["column_name", "field_name", "name"] if c in schema_df.columns), None)
if name_col is None:
    raise ValueError(f"Could not infer schema column-name field from DESCRIBE output columns: {schema_df.columns}")

input_columns = schema_df[name_col].to_list()
passthrough_cols = [c for c in input_columns if c not in ["genotypes"]]

missing_info = [f for f in ["AC", "AN"] if f not in input_columns]
if missing_info:
    print(
        f"Error: the allele frequency prior needs INFO/{' and INFO/'.join(missing_info)}, which "
        "the input header does not declare. Run `bcftools +fill-tags -t AC,AN` on the input first."
    )
    exit(1)


def quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'

# note on "start + 1" below: polars-bio for unknown reasons modifies the start column by subtracting 1 from it - this is probably a bug. For now we add 1 in the line below as a workaround.
select_passthrough_sql = ", ".join("i." + quote_ident(col) if col != "start" else "CAST(i." + quote_ident(col) + " + 1 AS INTEGER UNSIGNED) AS start" for col in passthrough_cols)
print(f"Selecting passthrough columns: {select_passthrough_sql}")

# p = max((AC + 1/2) / (AN + 1), --ds-prior-min-af). The half allele added to each side keeps p
# strictly inside (0, 1). AC is Number=A, so it arrives as a list.
ALLELE_FREQ_SQL = (
    f'GREATEST((CAST(COALESCE("AC"[1], 0) AS DOUBLE) + 0.5)'
    f' / (CAST(COALESCE("AN", 0) AS DOUBLE) + 1.0), {args.ds_prior_min_af!r})'
)

# Prior weights pi_g. --pl-type is fixed for the run, so the choice is baked into the SQL instead
# of branched on per genotype. prior_alt = 0 for a haploid genotype, which is what reduces the
# diploid dosage below to the haploid one.
if pl_type == "posterior":
    # PL already carries the caller's own prior; normalising it is all that is wanted.
    PRIOR_REF_SQL, PRIOR_HET_SQL = "1.0", "1.0"
    PRIOR_ALT_SQL = "CASE WHEN ploidy = 1 THEN 0.0 ELSE 1.0 END"
else:
    # Hardy-Weinberg genotype frequencies at allele frequency p.
    PRIOR_REF_SQL = (
        "CASE WHEN ploidy = 1 THEN 1.0 - allele_freq"
        " ELSE (1.0 - allele_freq) * (1.0 - allele_freq) END"
    )
    PRIOR_HET_SQL = (
        "CASE WHEN ploidy = 1 THEN allele_freq ELSE 2.0 * allele_freq * (1.0 - allele_freq) END"
    )
    PRIOR_ALT_SQL = "CASE WHEN ploidy = 1 THEN 0.0 ELSE allele_freq * allele_freq END"

DS_CTES_SQL = f"""\
WITH indexed AS (
  SELECT
    *,
    ROW_NUMBER() OVER () AS variant_idx,
    -- Repeated once per sample so it is an array column like GT: the only shape that keeps the
    -- FusedArrayTransform optimisation. GT is used for the count because it is the one FORMAT
    -- field every genotype carries.
    array_repeat({ALLELE_FREQ_SQL}, array_length(genotypes."GT")) AS allele_freq_arr
  FROM vcf_table
),
samples_unnested AS (
  SELECT
    variant_idx,
    UNNEST(genotypes."GT") AS gt,
    UNNEST(genotypes."GQ") AS gq,
    UNNEST(genotypes."DP") AS dp,
    UNNEST(genotypes."PL") AS pl,
    UNNEST(allele_freq_arr) AS allele_freq
  FROM indexed
),
sample_ploidy AS (
  SELECT
    *,
    -- Separators + 1. A missing haploid GT arrives as NULL, so ploidy is NULL for it.
    character_length(gt) - character_length(replace(replace(gt, '/', ''), '|', '')) + 1 AS ploidy
  FROM samples_unnested
),
sample_processed AS (
  SELECT
    *,
    (gq IS NOT NULL AND dp IS NOT NULL) AS is_good,
    (gt IN ('0/0', '0|0', '0')) AS is_hom_ref,
    COALESCE(pl[1], 0) AS pl0,
    COALESCE(pl[2], 0) AS pl1,
    COALESCE(pl[3], 0) AS pl2,
    -- Records are biallelic here (BCFTOOLS_NORM ran), so PL must carry ploidy+1 values.
    (pl IS NOT NULL AND array_length(pl) = ploidy + 1) AS pl_valid,
    (gq IS NOT NULL AND gq >= {args.calc_ds_min_gq}) AS gq_sufficient
  FROM sample_ploidy
),
reconstruction AS (
  SELECT
    *,
    -- 0 for a well-formed PL. Above 0 the sample's own genotype is not one of the ploidy+1 this
    -- record can express, which is what a multiallelic split leaves behind; see the docstring.
    CASE WHEN ploidy = 1 THEN LEAST(pl0, pl1) ELSE LEAST(pl0, pl1, pl2) END AS pl_min,
    -- Copies of this record's ALT according to GT. Spelled out rather than counted out of the
    -- string: a replace() over gt inside float arithmetic loses FusedArrayTransform (see the
    -- work plan). After the biallelic split these are all the genotypes there are.
    CASE gt
      WHEN '0/1' THEN 1.0 WHEN '1/0' THEN 1.0 WHEN '0|1' THEN 1.0 WHEN '1|0' THEN 1.0
      WHEN '1/1' THEN 2.0 WHEN '1|1' THEN 2.0 WHEN '1' THEN 1.0
      ELSE 0.0
    END AS gt_alt_count,
    -- An all-zero or absent PL says nothing about which genotype is right, so GQ is all there is.
    (pl IS NULL OR (pl_valid AND pl0 = 0 AND pl1 = 0 AND pl2 = 0)) AS pl_uninformative,
    -- w = P(the call is wrong), capped at ploidy/(ploidy+1): a called genotype is the most
    -- probable one, so it keeps at least 1/(ploidy+1).
    LEAST(POWER(10.0, -gq / 10.0), CASE WHEN ploidy = 1 THEN 0.5 ELSE 2.0 / 3.0 END)
      AS call_err_prob
  FROM sample_processed
),
with_ds AS (
  SELECT
    *,
    CASE WHEN pl0 < 255 THEN POWER(10.0, -pl0/10.0) ELSE 0.0 END AS likelihood_ref,
    CASE WHEN pl1 < 255 THEN POWER(10.0, -pl1/10.0) ELSE 0.0 END AS likelihood_het,
    CASE WHEN pl2 < 255 THEN POWER(10.0, -pl2/10.0) ELSE 0.0 END AS likelihood_alt,
    {PRIOR_REF_SQL} AS prior_ref,
    {PRIOR_HET_SQL} AS prior_het,
    {PRIOR_ALT_SQL} AS prior_alt
  FROM reconstruction
),
ds_terms AS (
  SELECT
    *,
    -- Unnormalised posteriors pi_g * L_g; the normalisation is the division in ds_value below.
    prior_ref * likelihood_ref AS post_ref,
    prior_het * likelihood_het AS post_het,
    prior_alt * likelihood_alt AS post_alt,
    -- P(the truth is one of the genotypes this record can express). 1 for a well-formed PL.
    POWER(10.0, -pl_min / 10.0) AS likelihood_weight
  FROM with_ds
),
ds_blend AS (
  SELECT
    *,
    -- DS = sum_g g * pi_g * L_g / sum_g pi_g * L_g, weighted by the mass this record accounts
    -- for; the rest belongs to genotypes it cannot express, whose copies of this ALT only GT
    -- knows. A well-formed PL has likelihood_weight = 1 and the second term vanishes. prior_alt
    -- is 0 for a haploid genotype, which reduces the first term to its haploid form.
    CASE WHEN (post_ref + post_het + post_alt) > 0 THEN
      likelihood_weight * (post_het + 2.0 * post_alt) / (post_ref + post_het + post_alt)
      + (1.0 - likelihood_weight) * gt_alt_count
    END AS ds_value
  FROM ds_terms
),
ds_calc AS (
  SELECT
    *,
    CASE
      WHEN gt IS NULL OR is_good THEN gt
      WHEN ploidy = 1 THEN '.'
      ELSE './.'
    END AS gt_final,
    CASE
      -- w split over the other genotypes by their Hardy-Weinberg frequencies, giving 2w/(2-p);
      -- see the module docstring. Only hom-ref calls qualify: for a carrier, GQ does not say
      -- which genotype the truth would be.
      WHEN pl_uninformative AND is_hom_ref AND gq_sufficient THEN
        CASE
          WHEN ploidy = 1 THEN call_err_prob
          ELSE 2.0 * call_err_prob / (2.0 - allele_freq)
        END
      WHEN pl_uninformative THEN NULL
      -- Too far from the genotypes this record can express to say anything: see the docstring.
      WHEN pl_min >= {args.ds_max_pl_min} THEN NULL
      WHEN pl_valid AND gq_sufficient AND ds_value IS NOT NULL THEN
        CASE WHEN ds_value < 0.000001 THEN 0.0 ELSE ds_value END
      ELSE NULL
    END AS ds
  FROM ds_blend
)"""

FILTER_SQL = DS_CTES_SQL + f""",
genotypes_aggregated AS (
  SELECT
    variant_idx,
    STRUCT(
      array_agg(gt_final) AS GT,
      array_agg(gq) AS GQ,
      array_agg(dp) AS DP,
      -- The rebuilt PL feeds DS only; the caller's own PL is what gets written.
      array_agg(pl) AS PL,
      array_agg(CAST(ds AS FLOAT)) AS DS
    ) AS genotypes
  FROM ds_calc
  GROUP BY variant_idx
)
SELECT {select_passthrough_sql},
  g.genotypes
FROM indexed i join genotypes_aggregated g ON i.variant_idx = g.variant_idx
"""

# The tracking counts, over the same CTEs. This is a second pass over the input rather than a
# by-product of the first: the query above has to end in the genotypes struct that sink_vcf
# writes, and nothing can be aggregated out of it on the way.
COUNTS_SQL = DS_CTES_SQL + """
SELECT
  COUNT(*) AS total,
  0 AS passed_through,
  COUNT(*) FILTER (WHERE ds IS NOT NULL AND NOT pl_uninformative) AS computed,
  COUNT(*) FILTER (WHERE ds IS NOT NULL AND pl_uninformative) AS reconstructed,
  COUNT(*) FILTER (WHERE ds IS NULL) AS missing
FROM ds_calc
"""


lf = pb.sql("EXPLAIN select * FROM (" + FILTER_SQL + ")")

optim_used = False
with pl.Config(fmt_str_lengths=None, tbl_rows=-1, tbl_cols=-1, fmt_table_cell_list_len=-1, tbl_width_chars=20000):
    logical_plan = lf.collect()["plan"][0]
    physical_plan = lf.collect()["plan"][1]
    print(f"sql result:")
    print("")
    print("  logical plan:")
    print(logical_plan)
    print("")
    print("  physical plan:")
    print(physical_plan)
    if "FusedArrayTransform" in physical_plan:
        optim_used = True

print("")
print(f"optim_used = {optim_used}")

if not optim_used:
    print("Error: for some reason FusedArrayTransform optimization is not used. This optimization is mandatory for this query. Interrupting..")
    exit(1)


# Keep only the FORMAT fields present in the output genotypes struct
output_format_fields = {}
for field_name in ["GT", "GQ", "DP", "PL"]:
    if field_name in input_format_fields:
        output_format_fields[field_name] = input_format_fields[field_name]

# Reuse the input's own DS definition when it has one, so DS is never declared twice.
output_format_fields["DS"] = input_ds or {
    # Number=1: BCFTOOLS_NORM split the record upstream, so there is exactly one ALT to dose.
    # The Number=A warning above is about input files, which have not been through that split.
    "number": "1",
    "type": "Float",
    "description": "Genotype dosage from PL fields",
}

header["format_fields"] = output_format_fields


lf = pb.sql(FILTER_SQL)
pb.set_source_metadata(lf, format="vcf", path=args.input_vcf_path, header=header)

pb.sink_vcf(lf, args.output_vcf_path)
print(f"Wrote: {args.output_vcf_path}")

write_tracking(passthrough=False)
