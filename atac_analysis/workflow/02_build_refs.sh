#!/bin/bash
set -euo pipefail

PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
REFDIR="$PROJECT/refs"
mkdir -p "$REFDIR"

# --- Build hg38  ---
cd "$REFDIR"

if [ ! -f hg38.fa ]; then
  echo "Downloading hg38 reference FASTA..."
  curl -L -o hg38.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
  gunzip hg38.fa.gz
fi

if [ ! -f hg38.fa.bwt ]; then
  echo "Building BWA index for hg38..."
  bwa index hg38.fa
fi

if [ ! -f hg38.chrom.sizes ]; then
  echo "Creating chrom.sizes..."
  faidx hg38.fa
  cut -f1,2 hg38.fa.fai > hg38.chrom.sizes
fi

# Optional: blacklist
if [ ! -f hg38.blacklist.bed ]; then
  echo "Downloading ENCODE blacklist (hg38)..."
  curl -L -o hg38.blacklist.bed https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
  gunzip hg38.blacklist.bed.gz
fi

echo " hg38 reference ready."
