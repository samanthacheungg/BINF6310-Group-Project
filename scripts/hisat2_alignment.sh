#!/usr/bin/env bash
set -euo pipefail

mkdir -p bam

# Go to project root or raw_data
cd "$(dirname "$0")/../raw_data"

for r1 in fastq/*_1.fastq; do
    base=$(basename "$r1" _1.fastq)
    r2="fastq/${base}_2.fastq"

    echo "Aligning $base ..."

    hisat2 \
      -x ref/yeast_ERCC_index \
      -1 "$r1" \
      -2 "$r2" \
      -S "bam/${base}.sam" \
      --phred33

    samtools view -bS "bam/${base}.sam" | \
      samtools sort -o "bam/${base}_sorted.bam" -

    samtools index "bam/${base}_sorted.bam"
    rm "bam/${base}.sam"
done

