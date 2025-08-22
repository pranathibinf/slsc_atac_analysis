#!/usr/bin/env bash
# Minimal ATAC-seq pipeline (post-refs) with robust resume behavior.
# Steps: fastp → bwa mem → name-sort → fixmate → coord-sort → markdup → MACS → FRiP
set -euo pipefail

PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
THREADS="${THREADS:-8}"
FORCE="${FORCE:-0}"   # 1 = redo even if outputs exist

REFDIR="$PROJECT/refs"
SUB="$PROJECT/data/processed"
WORK="$PROJECT/work"
ALIGN="$PROJECT/align"
PEAKS="$PROJECT/peaks"
QC="$PROJECT/qc"

SAMPLES=("SRX26680000" "SRX26680002" "SRX26680004")  # NeverMet, Met, FOXA2+ PDX

mkdir -p "$WORK" "$ALIGN" "$PEAKS" "$QC"

need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 1; }; }
need fastp; need bwa; need samtools; need bedtools
# MACS detect
MACS_BIN=""
if command -v macs3 >/dev/null 2>&1; then MACS_BIN="macs3";
elif command -v macs2 >/dev/null 2>&1; then MACS_BIN="macs2";
else echo "ERROR: need macs3 or macs2 installed"; exit 1; fi

genome_for() {
  case "$1" in
    SRX26680004) echo "mm10" ;;  # PDX mouse
    *)           echo "hg38" ;;  # FFPE human
  esac
}

# ------------------ 1) TRIM/QC ------------------
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

# ----------- 2) ALIGN → nsort → fixmate → sort → markdup -----------
echo "=== [2/4] bwa mem → fixmate chain ==="

clean_align_outputs () {
  local ACC="$1"
  rm -f "$ALIGN/${ACC}.nsort.bam" \
        "$ALIGN/${ACC}.fixmate.bam" \
        "$ALIGN/${ACC}.sorted.bam" "$ALIGN/${ACC}.sorted.bam.bai" \
        "$ALIGN/${ACC}.dedup.bam"  "$ALIGN/${ACC}.dedup.bam.bai" 2>/dev/null || true
}

for ACC in "${SAMPLES[@]}"; do
  REFBASE="$(genome_for "$ACC")"
  REFPREFIX="$REFDIR/${REFBASE}"  # bwa index prefix (no .fa)
  [[ -f "${REFPREFIX}.bwt" ]] || { echo "Missing BWA index for $REFBASE at ${REFPREFIX}.*"; exit 1; }

  TRIM1="$WORK/${ACC}_R1.trim.fastq.gz"
  TRIM2="$WORK/${ACC}_R2.trim.fastq.gz"

  NSORT="$ALIGN/${ACC}.nsort.bam"
  FIXM="$ALIGN/${ACC}.fixmate.bam"
  SORTB="$ALIGN/${ACC}.sorted.bam"
  DEDUP="$ALIGN/${ACC}.dedup.bam"

  if [[ "$FORCE" -eq 1 ]]; then clean_align_outputs "$ACC"; fi
  if [[ -f "$DEDUP" && -f "$DEDUP.bai" ]]; then
    echo "$ACC: dedup BAM exists -> skip alignment"
    continue
  fi

  if [[ ! -f "$NSORT" ]]; then
    echo "$ACC: bwa mem → name-sort"
    bwa mem -t "$THREADS" "$REFPREFIX" "$TRIM1" "$TRIM2" \
      | samtools view -b - \
      | samtools sort -n -@ "$THREADS" -o "$NSORT" -
  else
    echo "$ACC: using existing name-sorted BAM"
  fi

  if [[ ! -f "$FIXM" ]]; then
    echo "$ACC: samtools fixmate"
    samtools fixmate -m "$NSORT" "$FIXM"
  else
    echo "$ACC: using existing fixmate BAM"
  fi

  echo "$ACC: coordinate sort (rebuild from fixmate)"
  samtools sort -@ "$THREADS" -o "$SORTB" "$FIXM"
  samtools index "$SORTB"

  echo "$ACC: markdup (remove PCR dups)"
  samtools markdup -r "$SORTB" "$DEDUP"
  samtools index "$DEDUP"

  rm -f "$NSORT" "$FIXM"
done

# ------------------ 3) PEAKS + blacklist ------------------
echo "=== [3/4] $MACS_BIN peaks + blacklist ==="
call_and_filter () {
  local ACC="$1" GSIZE="$2" BLFILE="$3"
  local RAW="$PEAKS/${ACC}_peaks.narrowPeak"
  local FILT="$PEAKS/${ACC}.peaks.filt.bed"
  local BAM="$ALIGN/${ACC}.dedup.bam"

  if [[ "$FORCE" -eq 1 || ! -f "$RAW" ]]; then
    echo "$ACC: calling peaks"
    "$MACS_BIN" callpeak -t "$BAM" \
      -f BAMPE -g "$GSIZE" -n "$ACC" --outdir "$PEAKS" \
      --nomodel --shift -100 --extsize 200 --keep-dup all --pvalue 1e-3 >/dev/null
  else
    echo "$ACC: peaks exist -> skip callpeak"
  fi

  if [[ "$FORCE" -eq 1 || ! -f "$FILT" ]]; then
    if [[ -f "$REFDIR/$BLFILE" ]]; then
      echo "$ACC: blacklist filter ($BLFILE)"
      bedtools intersect -v -a "$RAW" -b "$REFDIR/$BLFILE" > "$FILT"
    else
      echo "WARN: $BLFILE missing; copying raw peaks."
      cp -f "$RAW" "$FILT"
    fi
  else
    echo "$ACC: filtered peaks exist -> skip"
  fi
}
call_and_filter "SRX26680000" hs "hg38.blacklist.bed"
call_and_filter "SRX26680002" hs "hg38.blacklist.bed"
call_and_filter "SRX26680004" mm "mm10.blacklist.bed"

# ------------------ 4) FRiP (robust) ------------------
echo "=== [4/4] FRiP ==="
OUTTSV="$QC/frip.tsv"
echo -e "Sample\tInPeaks\tTotal\tFRiP" > "$OUTTSV"
for ACC in "${SAMPLES[@]}"; do
  BAM="$ALIGN/${ACC}.dedup.bam"
  FRAGS="$ALIGN/${ACC}.frags.bed"

  # proper pair, primary, mapped, not QC-fail, not dup
  samtools view -@ "$THREADS" -b -f 0x2 -F 0x904 "$BAM" \
  | bedtools bamtobed -bedpe -i - \
  | awk 'BEGIN{OFS="\t"} ($1==$4) {s=($2<$5)?$2:$5; e=($3>$6)?$3:$6; if (s>=0 && e>s) print $1,s,e; }' \
  > "$FRAGS"

  TOTAL=$(wc -l < "$FRAGS")
  OVER=$(bedtools intersect -u -a "$FRAGS" -b "$PEAKS/${ACC}.peaks.filt.bed" | wc -l)
  awk -v a="$ACC" -v o="$OVER" -v t="$TOTAL" 'BEGIN{printf "%s\t%d\t%d\t%.4f\n", a,o,t,(t?o/t:0)}' >> "$OUTTSV"
done

echo "FRiP -> $OUTTSV"
echo "Done ✅"
