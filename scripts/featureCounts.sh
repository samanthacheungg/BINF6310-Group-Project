#!/usr/bin/env bash
set -euo pipefail

# Always run from the directory where this script lives (project root)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Paths (relative to project root)
ANNOTATION="raw_data/ref/new_yeast_R64_ERCC_refnames.gtf"
COUNTS_DIR="counts"
BAM_DIR="alignment/bam"

# Make sure output directory exists
mkdir -p "${COUNTS_DIR}"

echo "Using annotation: ${ANNOTATION}"
echo "Using BAM files in: ${BAM_DIR}"
echo "Output will be written to: ${COUNTS_DIR}/gene_counts.txt"
echo

# Run featureCounts
featureCounts \
  -a "${ANNOTATION}" \
  -o "${COUNTS_DIR}/gene_counts.txt" \
  -T 4 \
  -g gene_id \
  -t exon \
  -p \
  "${BAM_DIR}"/*_sorted.bam

echo
echo "featureCounts completed successfully."
echo "Counts file: ${COUNTS_DIR}/gene_counts.txt"

