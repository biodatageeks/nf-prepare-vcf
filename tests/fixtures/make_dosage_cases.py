#!/usr/bin/env python3
"""Generate the crafted dosage-edge-case fixtures for CALC_DOSAGE_POLARSBIO.

The big fixture (unprepared_rand_500.vcf.gz) is real 1000G data and is good at proving the
happy path, but it contains none of the cases the dosage code actually branches on: it has no
all-zero PL, no missing PL, no pre-existing DS, and its haploid genotypes were silently dropped
by the old implementation without any test noticing. Its multiallelic sites do produce split
records whose PL has no zero, but only at a handful of loci and never at a chosen min(PL).

This generator writes a small, hand-specified VCF where every genotype is one named case, so a
test can assert per-genotype behaviour by sample name instead of by aggregate statistics.

Usage:
    python3 make_dosage_cases.py OUTDIR

Writes five uncompressed VCFs into OUTDIR. All share one body; only the header differs:

    dosage_cases_no_ds.vcf     no DS declared, no DS values   -> pipeline must compute DS
    dosage_cases_ds_a.vcf      DS declared Number=A (correct) -> pipeline must pass DS through
    dosage_cases_ds_1.vcf      DS declared Number=1 (wrong)   -> pass through, but warn
    dosage_cases_caller_gatk.vcf          GATK FORMAT fingerprint -> --pl-type auto picks HWE
    dosage_cases_caller_deepvariant.vcf   DeepVariant fingerprint -> --pl-type auto picks flat

Compress/index them with bgzip + bcftools index -t.
"""

import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Samples. 24 is enough that AC/AN give a non-degenerate allele-frequency
# estimate, small enough that every genotype below can be written out by hand.
# Sex matters only on chrX outside the PAR: M* are haploid there, F* diploid.
# ---------------------------------------------------------------------------
MALES = [f"M{i:02d}" for i in range(1, 13)]
FEMALES = [f"F{i:02d}" for i in range(1, 13)]
SAMPLES = MALES + FEMALES

HEADER_COMMON = """\
##fileformat=VCFv4.2
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
"""

DS_HEADER = {
    "no_ds": "",
    "ds_a": '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Genotype dosage">\n',
    "ds_1": '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Genotype dosage">\n',
}

# Caller fingerprints for the `--pl-type auto` detection. These declare FORMAT fields that only
# one of the two callers emits; no genotype carries them, because detection reads the header
# alone. A caller-named ##source line is not used here: the headers that need fingerprinting in
# practice are exactly the ones that have no such line.
CALLER_HEADER = {
    "gatk": (
        '##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence">\n'
        '##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information">\n'
    ),
    "deepvariant": (
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">\n'
        '##FORMAT=<ID=MED_DP,Number=1,Type=Integer,Description="Median DP observed within the GVCF block">\n'
    ),
}

# ---------------------------------------------------------------------------
# Genotype cases.
#
# Each entry is (GT, DP, GQ, PL, DS). Use None for a missing value; it is
# written as '.'. DS is only emitted into the ds_a / ds_1 fixtures.
# ---------------------------------------------------------------------------

# --- confident, unambiguous diploid calls: DS must land on the integer ---
CONF_HOM_REF = ("0/0", 30, 99, [0, 120, 1200], 0.0)
CONF_HET = ("0/1", 28, 99, [900, 0, 950], 1.0)
CONF_HOM_ALT = ("1/1", 26, 99, [1400, 130, 0], 2.0)

# --- intermediate diploid calls: DS must be strictly between integers ---
WEAK_HET = ("0/1", 6, 20, [20, 0, 130], 1.0)
WEAK_HOM_REF = ("0/0", 5, 15, [0, 15, 180], 0.03)

# --- all-zero PL on a hom-ref call, spanning the --calc-ds-min-gq=3 threshold ---
AZ_GQ0 = ("0/0", 8, 0, [0, 0, 0], None)  # below threshold -> DS missing
AZ_GQ2 = ("0/0", 8, 2, [0, 0, 0], None)  # below threshold -> DS missing
AZ_GQ3 = ("0/0", 8, 3, [0, 0, 0], None)  # at threshold    -> DS reconstructed
AZ_GQ10 = ("0/0", 9, 10, [0, 0, 0], None)
AZ_GQ40 = ("0/0", 31, 40, [0, 0, 0], None)

# --- PL absent entirely, GQ present: treated like all-zero PL (decision 13.5) ---
NOPL_GQ20 = ("0/0", 12, 20, None, None)
NOPL_NOGQ = ("0/0", 12, None, None, None)  # nothing to go on -> DS missing, GT masked

# --- masking: GQ or DP missing means the call itself is not trustworthy ---
NO_GQ = ("0/1", 14, None, [40, 0, 200], None)
NO_DP = ("0/1", None, 30, [40, 0, 200], None)

# --- fully missing ---
MISSING_DIP = ("./.", None, None, None, None)

# --- pre-existing DS values, exercised only in the ds_a / ds_1 fixtures ---
# An imputed-looking dosage that disagrees with PL: proves DS is passed through
# untouched rather than recomputed.
DS_DISAGREES = ("0/0", 30, 99, [0, 120, 1200], 1.75)
# Out of range. Decision 13.7: if DS is present we trust it, so this is passed
# through as-is too -- the test pins that, so a future change to clamp it fails loudly.
DS_OUT_OF_RANGE = ("0/0", 30, 99, [0, 120, 1200], 2.5)
# DS missing on a genotype that does have usable PL: must be computed, not left missing.
DS_GAP = ("0/1", 28, 99, [900, 0, 950], None)

# --- haploid (chrX outside the PAR, male samples): 2 PL values, DS in [0,1] ---
HAP_CONF_REF = ("0", 30, 99, [0, 700], 0.0)
HAP_CONF_ALT = ("1", 28, 99, [800, 0], 1.0)
HAP_WEAK_ALT = ("1", 5, 12, [12, 0], 0.94)
HAP_AZ_GQ10 = ("0", 7, 10, [0, 0], None)  # haploid all-zero PL
HAP_AZ_GQ2 = ("0", 7, 2, [0, 0], None)  # below threshold
HAP_MISSING = (".", None, None, None, None)  # must stay '.', not become './.'

# --- PL length inconsistent with the ploidy taken from GT: DS must stay missing ---
BAD_PL_DIP = ("0/1", 20, 40, [30, 0], None)  # diploid GT, 2 PLs
BAD_PL_HAP = ("0", 20, 40, [0, 40, 400], None)  # haploid GT, 3 PLs

# --- PL with no zero in it, as a multiallelic split leaves behind ---------------------
# The sample's own genotype uses an allele this record does not carry, so the entry that was 0
# belongs to a genotype the record cannot express. What is left says which of the three
# representable genotypes is least bad, which is not the same question. min(PL) measures how far
# outside the record the truth sits: at or above --ds-max-pl-min the genotype gets no dosage,
# below it the dosage is weighted between the likelihoods and GT.
SPLIT_FAR_NONCARRIER = ("0/0", 25, 40, [331, 260, 247], None)  # min 247, hom-alt for another ALT
SPLIT_FAR_CARRIER = ("0/1", 24, 35, [300, 90, 60], None)  # min 60, compound het
SPLIT_AT_LIMIT = ("0/0", 20, 25, [50, 35, 30], None)  # min 30: at the default, dropped
SPLIT_BELOW_LIMIT = ("0/0", 20, 25, [49, 34, 29], None)  # min 29: just inside, dosed
SPLIT_BAND_NONCARRIER = ("0/0", 22, 30, [40, 25, 20], None)  # min 20 -> DS near 0
SPLIT_BAND_CARRIER = ("0/1", 21, 28, [45, 30, 25], None)  # min 25 -> DS near 1
# Haploid: the smallest of two values, not three. Taking the minimum over a diploid-shaped PL
# would always find the padding 0 and quietly exempt every haploid genotype from all of this.
SPLIT_HAP_FAR = ("0", 20, 30, [60, 40], None)  # min 40, dropped
SPLIT_HAP_BAND_ALT = ("1", 20, 30, [35, 20], None)  # min 20 -> DS near 1
SPLIT_HAP_BAND_REF = ("0", 20, 30, [20, 35], None)  # min 20 -> DS near 0


def rec(chrom, pos, vid, ref, alt, per_sample):
    """per_sample: dict sample -> case tuple. Samples not listed get a default."""
    return dict(chrom=chrom, pos=pos, vid=vid, ref=ref, alt=alt, per_sample=per_sample)


def build_records():
    records = []

    # -- r1: common autosomal variant. Confident calls of all three genotypes, so
    #    AC/AN are substantial and the HWE prior is close to flat. DS must equal the
    #    GT allele count for every sample here.
    ps = {}
    for s in SAMPLES[:8]:
        ps[s] = CONF_HOM_REF
    for s in SAMPLES[8:18]:
        ps[s] = CONF_HET
    for s in SAMPLES[18:]:
        ps[s] = CONF_HOM_ALT
    records.append(rec("chr22", 1000, "common", "A", "G", ps))

    # -- r2: singleton. One weak het carrier against 23 confident hom-refs. This is the
    #    record where the flat and HWE priors disagree most (analysis section 7.3): the
    #    flat prior hands ~0.97 dosage to the carrier and spurious dosage to non-carriers.
    ps = {s: CONF_HOM_REF for s in SAMPLES}
    ps[SAMPLES[0]] = WEAK_HET
    ps[SAMPLES[1]] = WEAK_HOM_REF
    ps[SAMPLES[2]] = WEAK_HOM_REF
    records.append(rec("chr22", 2000, "singleton", "C", "T", ps))

    # -- r3: all-zero PL across the GQ threshold.
    ps = {s: CONF_HOM_REF for s in SAMPLES}
    ps[SAMPLES[0]] = AZ_GQ0
    ps[SAMPLES[1]] = AZ_GQ2
    ps[SAMPLES[2]] = AZ_GQ3
    ps[SAMPLES[3]] = AZ_GQ10
    ps[SAMPLES[4]] = AZ_GQ40
    ps[SAMPLES[5]] = CONF_HET  # keep AC > 0 so p-hat is not at its floor
    ps[SAMPLES[6]] = CONF_HET
    records.append(rec("chr22", 3000, "allzero_pl", "G", "A", ps))

    # -- r4: missing PL / missing GQ / missing DP / fully missing.
    ps = {s: CONF_HOM_REF for s in SAMPLES}
    ps[SAMPLES[0]] = NOPL_GQ20
    ps[SAMPLES[1]] = NOPL_NOGQ
    ps[SAMPLES[2]] = NO_GQ
    ps[SAMPLES[3]] = NO_DP
    ps[SAMPLES[4]] = MISSING_DIP
    ps[SAMPLES[5]] = CONF_HET
    records.append(rec("chr22", 4000, "missing_fields", "T", "C", ps))

    # -- r5: pre-existing DS, including a value that disagrees with PL and one out of range.
    ps = {s: CONF_HOM_REF for s in SAMPLES}
    ps[SAMPLES[0]] = DS_DISAGREES
    ps[SAMPLES[1]] = DS_OUT_OF_RANGE
    ps[SAMPLES[2]] = DS_GAP
    ps[SAMPLES[3]] = CONF_HET
    records.append(rec("chr22", 5000, "existing_ds", "A", "T", ps))

    # -- r6: PL length inconsistent with GT ploidy.
    ps = {s: CONF_HOM_REF for s in SAMPLES}
    ps[SAMPLES[0]] = BAD_PL_DIP
    records.append(rec("chr22", 6000, "bad_pl_len", "C", "G", ps))

    # -- r7: what one record of a split multiallelic site looks like -- genotypes whose PL has
    #    no zero, spanning the --ds-max-pl-min threshold on both sides.
    ps = {s: CONF_HOM_REF for s in SAMPLES}
    ps[SAMPLES[0]] = SPLIT_FAR_NONCARRIER
    ps[SAMPLES[1]] = SPLIT_FAR_CARRIER
    ps[SAMPLES[2]] = SPLIT_AT_LIMIT
    ps[SAMPLES[3]] = SPLIT_BELOW_LIMIT
    ps[SAMPLES[4]] = SPLIT_BAND_NONCARRIER
    ps[SAMPLES[5]] = SPLIT_BAND_CARRIER
    ps[SAMPLES[6]] = CONF_HET  # keeps AC > 0, so p-hat is not at its floor
    ps[SAMPLES[7]] = CONF_HET
    records.append(rec("chr22", 7000, "split_no_zero_pl", "A", "T", ps))

    # -- r8: chrX outside the PAR. Males haploid, females diploid, at one record.
    #    AN must count 1 allele per male and 2 per female.
    ps = {}
    for s in MALES[:5]:
        ps[s] = HAP_CONF_REF
    for s in MALES[5:8]:
        ps[s] = HAP_CONF_ALT
    ps[MALES[8]] = HAP_WEAK_ALT
    ps[MALES[9]] = HAP_AZ_GQ10
    ps[MALES[10]] = HAP_AZ_GQ2
    ps[MALES[11]] = HAP_MISSING
    for s in FEMALES[:6]:
        ps[s] = CONF_HOM_REF
    for s in FEMALES[6:10]:
        ps[s] = CONF_HET
    for s in FEMALES[10:]:
        ps[s] = CONF_HOM_ALT
    records.append(rec("chrX", 7000, "x_nonpar_mixed", "A", "C", ps))

    # -- r9: haploid PL length mismatch, on its own record so it cannot perturb r8's AC/AN.
    ps = {s: HAP_CONF_REF for s in MALES}
    ps.update({s: CONF_HOM_REF for s in FEMALES})
    ps[MALES[0]] = BAD_PL_HAP
    records.append(rec("chrX", 8000, "x_bad_pl_len", "G", "T", ps))

    # -- r10: the same as r7 on chrX, where the males are haploid.
    ps = {s: HAP_CONF_REF for s in MALES}
    ps.update({s: CONF_HOM_REF for s in FEMALES})
    ps[MALES[0]] = SPLIT_HAP_FAR
    ps[MALES[1]] = SPLIT_HAP_BAND_ALT
    ps[MALES[2]] = SPLIT_HAP_BAND_REF
    ps[MALES[3]] = HAP_CONF_ALT
    ps[FEMALES[0]] = CONF_HET
    records.append(rec("chrX", 9000, "x_split_no_zero_pl", "C", "A", ps))

    return records


def fmt_cell(case, with_ds):
    gt, dp, gq, pl, ds = case
    parts = [
        gt,
        "." if dp is None else str(dp),
        "." if gq is None else str(gq),
        "." if pl is None else ",".join(str(v) for v in pl),
    ]
    if with_ds:
        parts.append("." if ds is None else f"{ds:g}")
    return ":".join(parts)


def ac_an(record):
    """AC/AN over the GT fields, exactly as `bcftools +fill-tags -t AC,AN` computes them.

    AN counts one allele per haploid and two per diploid non-missing genotype, so a chrX
    record shared by male and female samples contributes a mix of both.
    """
    ac = an = 0
    for s in SAMPLES:
        gt = record["per_sample"][s][0]
        for allele in gt.replace("|", "/").split("/"):
            if allele == ".":
                continue
            an += 1
            if allele != "0":
                ac += 1
    return ac, an


def write_vcf(path, ds_mode, caller_header=""):
    with_ds = ds_mode != "no_ds"
    fmt_keys = "GT:DP:GQ:PL" + (":DS" if with_ds else "")

    lines = (HEADER_COMMON + DS_HEADER[ds_mode] + caller_header).splitlines()
    lines.append("#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + SAMPLES))

    for r in build_records():
        cells = [fmt_cell(r["per_sample"][s], with_ds) for s in SAMPLES]
        ac, an = ac_an(r)
        info = f"AC={ac};AN={an}"
        lines.append(
            "\t".join([r["chrom"], str(r["pos"]), r["vid"], r["ref"], r["alt"], "100", "PASS", info, fmt_keys] + cells)
        )

    Path(path).write_text("\n".join(lines) + "\n")
    print(f"wrote {path}")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    outdir = Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)
    for mode in DS_HEADER:
        write_vcf(outdir / f"dosage_cases_{mode}.vcf", mode)
    for caller, caller_header in CALLER_HEADER.items():
        write_vcf(outdir / f"dosage_cases_caller_{caller}.vcf", "no_ds", caller_header)


if __name__ == "__main__":
    main()
