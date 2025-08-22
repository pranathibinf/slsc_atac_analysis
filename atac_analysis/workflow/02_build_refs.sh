#!/usr/bin/env bash
# Build BWA indexes for hg38 and mm10 and fetch ENCODE blacklists
set -euo pipefail

PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
REFS="$PROJECT/refs"
mkdir -p "$REFS"
cd "$REFS"

# ---- hg38 ----
if [ ! -f hg38.fa ]; then
  curl -L http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > hg38.fa
fi
bwa index -p hg38 hg38.fa
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.chrom.sizes

# ---- mm10 ----
if [ ! -f mm10.fa ]; then
  curl -L http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz | gunzip -c > mm10.fa
fi
bwa index -p mm10 mm10.fa
samtools faidx mm10.fa
cut -f1,2 mm10.fa.fai > mm10.chrom.sizes

# ---- ENCODE blacklists (v2) ----
curl -L -o hg38.blacklist.bed.gz https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
gunzip -f hg38.blacklist.bed.gz

curl -L -o mm10.blacklist.bed.gz https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz
gunzip -f mm10.blacklist.bed.gz

echo "Refs ready in: $REFS"
