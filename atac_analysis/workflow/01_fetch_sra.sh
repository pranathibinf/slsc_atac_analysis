#!/usr/bin/env bash
# Download FASTQs from SRA into data/raw/ and (optionally) write small subsamples to data/processed/.
# Requires: fastq-dl OR sra-toolkit (fasterq-dump), pigz, seqtk (optional for downsampling)
set -euo pipefail

PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
RAW="$PROJECT/data/raw"
PROC="$PROJECT/data/processed"

mkdir -p "$RAW" "$PROC"

# ---- Choose one method ----
# Method 1: fastq-dl (recommended if installed)
dl_acc() {
  local acc="$1"
  fastq-dl -a "$acc" -o "$RAW/${acc}"
}

# Method 2: fasterq-dump (uncomment to use)
# dl_acc() {
#   local acc="$1"
#   mkdir -p "$RAW/${acc}"
#   fasterq-dump --split-files --threads 8 -O "$RAW/${acc}" "$acc"
#   pigz -f "$RAW/${acc}"/*.fastq
# }

# SRX accessions (will expand to SRRs automatically with fastq-dl)
for ACC in SRX26680000 SRX26680002 SRX26680004; do
  echo "Downloading $ACC ..."
  dl_acc "$ACC"
done

echo "Raw FASTQs in: $RAW"

# Optional: create small subsamples (1M pairs) to data/processed/
# Requires seqtk
if command -v seqtk >/dev/null 2>&1; then
  echo "Creating 1M-pair subsamples in $PROC ..."
  for R1 in "$RAW"/*/*_1.fastq.gz "$RAW"/*/*R1*.fastq.gz 2>/dev/null; do
    [[ -f "$R1" ]] || continue
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    [[ -f "$R2" ]] || R2="${R1/R1/R2}"
    base=$(basename "$R1" | sed 's/_1\.fastq\.gz//; s/_R1\.fastq\.gz//')
    seqtk sample -s100 "$R1" 1000000 | pigz > "$PROC/${base}_R1.sub.fastq.gz"
    seqtk sample -s100 "$R2" 1000000 | pigz > "$PROC/${base}_R2.sub.fastq.gz"
  done
  echo "Processed subsamples in: $PROC"
else
  echo "NOTE: seqtk not found; skipping subsampling. Provide subsamples in $PROC manually."
fi
