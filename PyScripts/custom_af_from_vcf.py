#!/usr/bin/env python3

from cyvcf2 import VCF
import csv
import argparse

def aggregate_variant_data(vcf):
    variant_data = {}
    for variant in vcf:
        if "<*>" in variant.ALT:  # skip non-variant records
            continue
        position = variant.POS
        ref = variant.REF
        alts = variant.ALT
        # init data structure
        if position not in variant_data:
            variant_data[position] = []
        ad = variant.format("AD") if "AD" in variant.FORMAT else None
        variant_data[position].append((ref, alts, ad))
    return variant_data

def calculate_adjusted_AF(variant_data, samples, writer):
    for position, variants in variant_data.items():
        for ref, alts, ad in variants:
            row = {
                "Position": position,
                "Ref": ref,
                "Alt": ",".join(alts),
            }
            for i, sample_name in enumerate(samples):
                # init counts
                ref_count = alt_count = total_count = 0
                if ad is not None:
                    ref_count = ad[i][0]
                    alt_counts = ad[i][1:]
                    for alt_ad in alt_counts:
                        alt_count += alt_ad  # sum alt allele counts
                    total_count = ref_count + alt_count
                af = alt_count / total_count if total_count > 0 else 0
                row[sample_name] = "{:.4f}".format(af)
            writer.writerow(row)

def handle_gvcf_and_calculate_AF(vcf_file, output_csv):
    vcf = VCF(vcf_file)
    samples = vcf.samples
    fieldnames = ["Position", "Ref", "Alt"] + samples
    with open(output_csv, mode="w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        # aggregate variant data
        variant_data = aggregate_variant_data(vcf)
        # calulate AF considering other variants at the same position
        calculate_adjusted_AF(variant_data, samples, writer)

def main():
    parser = argparse.ArgumentParser(
        description="Calculate adjusted AF for variants in a VCF."
    )
    parser.add_argument("vcf_file", help="Path to the input VCF file (gzipped).")
    parser.add_argument("output_csv", help="Path to the output CSV file.")
    args = parser.parse_args()
    handle_gvcf_and_calculate_AF(args.vcf_file, args.output_csv)

if __name__ == "__main__":
    main()