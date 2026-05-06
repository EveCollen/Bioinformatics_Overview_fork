#!/bin/bash

# Exit on error
set -e

# Check at least one VCF provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <vcf1> [vcf2 ... vcfN]"
    exit 1
fi

# Store all input VCFs
vcfs=("$@")

echo "Running the validator script on each VCF..."
for vcf in "${vcfs[@]}"; do
    echo "Validating $vcf"
    python3 ./vcf_validator.py --input_vcf "$vcf"
done

echo "Running the awk script to pull out info about variants of interest..."
awk -f ./extract_variants.awk ../4_refs/DPYD_variants_genome_location.csv "${vcfs[@]}" > ../3_reports/variant_info.txt

echo "Generating report for QC checking..."
python3 ./generate_DPYD_report.py --input_vcfs "${vcfs[@]}" --ouput ../3_reports/Report_for_DPYD_variants.txt

echo "Done."
