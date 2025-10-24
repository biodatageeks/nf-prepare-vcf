#!/usr/bin/env python3

import cyvcf2
import pandas as pd
import sys
import warnings
import time
import os
import argparse

def current_milli_time():
    return round(time.time() * 1000)

def proccess_vcf(vcf_path, output_prefix, threads):
    gene_pos = {}  # gene -> (chrom, min_pos, max_pos)

    try:
        vcf = cyvcf2.VCF(vcf_path, threads=threads)
        for variant in vcf:
            csq = variant.INFO.get("CSQ")
            if not csq:
                continue
            csq_list = csq.split(",")
            seen_genes = set()  # Per variant to avoid duplicates
            for csq_entry in csq_list:
                fields = csq_entry.split("|")
                if len(fields) < 4:
                    continue
                gene = fields[3]
                if gene and gene not in seen_genes:
                    seen_genes.add(gene)
                    chrom = variant.CHROM
                    pos = variant.POS
                    if gene in gene_pos:
                        g_chrom, min_p, max_p = gene_pos[gene]
                        if g_chrom != chrom:
                            warnings.warn(f"Gene {gene} found on multiple chromosomes: {g_chrom} and {chrom}. Skipping.")
                            del gene_pos[gene]
                            continue
                        gene_pos[gene] = (chrom, min(min_p, pos), max(max_p, pos))
                    else:
                        gene_pos[gene] = (chrom, pos, pos)
    except Exception as e:
        sys.exit(f"Error parsing VCF: {e}")
    finally:
        vcf.close()

    # Create DataFrame
    data = [{'gene': gene, 'chrom': chrom, 'start': min_pos, 'end': max_pos} for gene, (chrom, min_pos, max_pos) in gene_pos.items()]
    df = pd.DataFrame(data)
    df.to_csv(output_prefix + '_gene_ranges.csv', index=False)

def main():
    parser = argparse.ArgumentParser(description="Process VCF to compute gene ranges.")
    parser.add_argument("--input-vcf", required=True, help="VCF.gz file (e.g., data.vcf.gz)")
    parser.add_argument("--output-prefix", help="Prefix for the output file")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads to use (default: 1)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_vcf):
        sys.exit(f"Error: Input VCF file {args.input_vcf} does not exist.")
    
    start_time = current_milli_time()
    # Process the VCF
    proccess_vcf(
        args.input_vcf,
        args.output_prefix,
        args.threads
    )
    print(f"Executed in {current_milli_time() - start_time} ms", flush=True)

if __name__ == "__main__":
    main()
