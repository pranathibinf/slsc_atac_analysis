#!/usr/bin/env bash
# Minimal ATAC-seq pipeline (post-refs) with robust resume behavior.
# Steps: fastp → bwa mem → name-sort → fixmate → coord-sort → markdup → MACS → FRiP
set -euo pipefail

# -------------------------- config --------------------------------
PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
THREADS="${THREADS:-8}"
FORCE="${FORCE:-0}"   # 1 = redo even if outputs exist

REFDIR="$PROJECT/refs"
SUB="$PROJECT/data/processed"
WORK="$PROJECT/work"
ALIGN="$PROJECT/align"
PEAKS="$PROJECT/peaks"
QC="$PROJECT/qc"

# Default samples 
SAMPLES=(SRX26680000)
# If SAMPLES is provided as a space-separated env var, convert to array
if [[ -n "${SAMPLES:-}" && "$SAMPLES" != *"("*")"* ]]; then
  read -r -a SAMPLES <<< "$SAMPLES"
fi

mkdir -p "$WORK" "$ALIGN" "$PEAKS" "$QC"

# -------------------------- tool checks ---------------------------------------
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 1; }; }
need fastp; need bwa; need samtools; need bedtools

# macs3 or macs2
MACS_BIN=""
if command -v macs3 >/dev/null 2>&1; then MACS_BIN="macs3";
elif command -v macs2 >/dev/null 2>&1; then MACS_BIN="macs2";
else echo "ERROR: need macs3 or macs2 installed"; exit 1; fi

# Map sample → genome 
genome_for() { echo "hg38"; }

# ------------------------------- helpers --------------------------------------
call_and_filter() {
  local ACC="$1"
  local genome_size_flag="$2"    # hs (hg38) / mm (mm10)
  local blacklist="$3"

  local BAM="$ALIGN/${ACC}.dedup.bam"
  local OUTPFX="$PEAKS/${ACC}"
  local RAWNP="${OUTPFX}_peaks.narrowPeak"
  local FILT="${PEAKS}/${ACC}.peaks.filt.bed"

  if [[ ! -f "$BAM" ]]; then
    echo "$ACC: skip (no BAM)"
    return 0
  fi

  if [[ "$FORCE" -eq 1 || ! -f "$RAWNP" ]]; then
    echo "$ACC: calling peaks"
    $MACS_BIN callpeak -t "$BAM" -f BAMPE -g "$genome_size_flag" -n "$ACC" --outdir "$PEAKS" \
      --pvalue 1e-3 --nomodel --shift -100 --extsize 200 --keep-dup all
  else
    echo "$ACC: peaks exist -> skip callpeak"
  fi

  if [[ -f "$RAWNP" && -f "$blacklist" ]]; then
    if [[ "$FORCE" -eq 1 || ! -f "$FILT" ]]; then
      bedtools intersect -v -a "$RAWNP" -b "$blacklist" > "$FILT"
      echo "$ACC: wrote blacklist-filtered peaks -> $FILT"
    else
      echo "$ACC: filtered peaks exist -> skip"
    fi
  fi
}

# === [1/4] fastp (trim/QC) ===
echo "=== [1/4] fastp (trim/QC) ==="
for ACC in "${SAMPLES[@]}"; do
  IN1="$SUB/${ACC}_R1.sub.fastq.gz"
  IN2="$SUB/${ACC}_R2.sub.fastq.gz"
  OUT1="$WORK/${ACC}_R1.trim.fastq.gz"
  OUT2="$WORK/${ACC}_R2.trim.fastq.gz"
  HTML="$QC/${ACC}.fastp.html"
  JSON="$QC/${ACC}.fastp.json"

  [[ -f "$IN1" && -f "$IN2" ]] || { echo "Missing $IN1 or $IN2"; exit 1; }

  if [[ "$FORCE" -eq 1 || ! -f "$OUT1" || ! -f "$OUT2" ]]; then
    echo "fastp -> $ACC"
    fastp -w "$THREADS" -i "$IN1" -I "$IN2" -o "$OUT1" -O "$OUT2" -h "$HTML" -j "$JSON" >/dev/null
  else
    echo "fastp -> $ACC (skip: outputs exist)"
  fi
done

# === [2/4] bwa mem → fixmate chain ===
echo "=== [2/4] bwa mem → fixmate chain ==="
for ACC in "${SAMPLES[@]}"; do
  REFBASE="$(genome_for "$ACC")"
  REFPREFIX="$REFDIR/${REFBASE}"  # IMPORTANT: prefix (no .fa)
  [[ -f "${REFPREFIX}.bwt" ]] || { echo "Missing BWA index for $REFBASE at ${REFPREFIX}.*"; exit 1; }

  TRIM1="$WORK/${ACC}_R1.trim.fastq.gz"
  TRIM2="$WORK/${ACC}_R2.trim.fastq.gz"
  NSORT="$ALIGN/${ACC}.nsort.bam"
  FIXM="$ALIGN/${ACC}.fixmate.bam"
  SORTB="$ALIGN/${ACC}.sorted.bam"
  DEDUP="$ALIGN/${ACC}.dedup.bam"

  if [[ "$FORCE" -eq 1 ]]; then
    rm -f "$NSORT" "$FIXM" "$SORTB" "$SORTB.bai" "$DEDUP" "$DEDUP.bai" 2>/dev/null || true
  fi

  if [[ -f "$DEDUP" && -f "$DEDUP.bai" ]]; then
    echo "$ACC: dedup BAM exists -> skip alignment"
    continue
  fi

  echo "$ACC: bwa mem"
  bwa mem -t "$THREADS" "$REFPREFIX" "$TRIM1" "$TRIM2" \
    | samtools view -bh -F 0x4 - \
    | samtools sort -n -@ "$THREADS" -o "$NSORT" -

  echo "$ACC: fixmate"
  samtools fixmate -m "$NSORT" "$FIXM"

  echo "$ACC: position sort"
  samtools sort -@ "$THREADS" -o "$SORTB" "$FIXM"
  samtools index "$SORTB"

  echo "$ACC: markdup (remove PCR dups)"
  samtools markdup -r -@ "$THREADS" "$SORTB" "$DEDUP"
  samtools index "$DEDUP"

  rm -f "$NSORT" "$FIXM" "$SORTB" "$SORTB.bai" 2>/dev/null || true
done

# === [3/4] macs3 peaks + blacklist ===
echo "=== [3/4] macs3 peaks + blacklist ==="
for ACC in "${SAMPLES[@]}"; do
  # hg38 only for now (PDX/mm10 can be added later if needed)
  call_and_filter "$ACC" hs "$REFDIR/hg38.blacklist.bed"
done

# === [4/4] FRiP ===
echo "=== [4/4] FRiP ==="
FRIP="$QC/frip.tsv"
echo -e "sample\treads_in_peaks\tproper_pairs\tfrip" > "$FRIP"

for ACC in "${SAMPLES[@]}"; do
  BAM="$ALIGN/${ACC}.dedup.bam"
  NP="${PEAKS}/${ACC}_peaks.narrowPeak"

  if [[ ! -f "$BAM" || ! -f "$NP" ]]; then
    echo "$ACC: skip FRiP (missing BAM or peaks)"
    continue
  fi

  PP=$(samtools view -c -f 0x2 -F 0x904 "$BAM")
  if [[ "$PP" -eq 0 ]]; then
    echo -e "${ACC}\t0\t0\t0.0" >> "$FRIP"
    continue
  fi

  QNAME="$ALIGN/${ACC}.qname.bam"
  FRAG="$ALIGN/${ACC}.frags.bed"
  samtools view -u -f 0x2 -F 0x904 "$BAM" \
    | samtools sort -n -@ "$THREADS" -o "$QNAME" -
  bedtools bamtobed -bedpe -i "$QNAME" \
    | awk '($1!="." && $2>=0 && $3>$2){print $1"\t"$2"\t"$3}' > "$FRAG"

  INPK=$(bedtools intersect -u -a "$FRAG" -b "$NP" | wc -l | tr -d ' ')
  FRIP_VAL=$(python3 - <<PY
pp=int("$PP"); ip=int("$INPK")
print(0.0 if pp==0 else round(ip/pp,4))
PY
)
  echo -e "${ACC}\t${INPK}\t${PP}\t${FRIP_VAL}" >> "$FRIP"
done

echo "FRiP -> $FRIP"
echo "Done ✅"
