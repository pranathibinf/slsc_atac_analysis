#!/usr/bin/env bash
set -euo pipefail

# -------------------- Config  --------------------
PROJECT="${PROJECT:-$(cd "$(dirname "$0")/.." && pwd)}"
THREADS="${THREADS:-8}"
FORCE="${FORCE:-0}"          # 1 = redo even if outputs exist
SAMPLES="${SAMPLES:-SRX26680000}"  # single sample id

# Layout
REFDIR="$PROJECT/refs"
SUB="$PROJECT/data/processed"
WORK="$PROJECT/work"
ALIGN="$PROJECT/align"
PEAKS="$PROJECT/peaks"
QC="$PROJECT/qc"

mkdir -p "$WORK" "$ALIGN" "$PEAKS" "$QC"

# -------------------- Tool checks --------------------
need(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 1; }; }
need fastp; need bwa; need samtools; need bedtools
MACS_BIN=""
if command -v macs3 >/dev/null 2>&1; then MACS_BIN="macs3";
elif command -v macs2 >/dev/null 2>&1; then MACS_BIN="macs2";
else echo "ERROR: need macs3 or macs2 installed"; exit 1; fi

# -------------------- Helpers --------------------
# Resolve BWA prefix: accept either refs/hg38.* or refs/hg38.fa.*
resolve_refprefix() {
  local base="$1"  # hg38
  if   [[ -f "$REFDIR/${base}.bwt" ]];     then echo "$REFDIR/${base}";
  elif [[ -f "$REFDIR/${base}.fa.bwt" ]];  then echo "$REFDIR/${base}.fa";
  else
    echo "ERROR: Cannot find BWA index for $base at $REFDIR/{${base}*.bwt}" >&2
    exit 1
  fi
}

# Human only for this repo (FFPE is human)
genome_for(){ echo "hg38"; }
gsize_for(){  echo "hs";   }  # MACS genome size token

# -------------------- 1) fastp trim/QC --------------------
echo "=== [1/4] fastp (trim/QC) ==="
for ACC in $SAMPLES; do
  IN1="$SUB/${ACC}_R1.sub.fastq.gz"
  IN2="$SUB/${ACC}_R2.sub.fastq.gz"
  OUT1="$WORK/${ACC}_R1.trim.fastq.gz"
  OUT2="$WORK/${ACC}_R2.trim.fastq.gz"
  HTML="$QC/${ACC}.fastp.html"
  JSON="$QC/${ACC}.fastp.json"

  [[ -s "$IN1" && -s "$IN2" ]] || { echo "Missing $IN1 or $IN2"; exit 1; }

  if [[ "$FORCE" -eq 1 || ! -s "$OUT1" || ! -s "$OUT2" || ! -s "$JSON" ]]; then
    echo "fastp -> $ACC"
    fastp -w "$THREADS" -i "$IN1" -I "$IN2" -o "$OUT1" -O "$OUT2" -h "$HTML" -j "$JSON" >/dev/null
  else
    echo "fastp -> $ACC (skip: outputs exist)"
  fi
done

# -------------------- 2) bwa mem → fixmate chain --------------------
echo "=== [2/4] bwa mem → fixmate chain ==="
for ACC in $SAMPLES; do
  REFBASE="$(genome_for "$ACC")"            # hg38
  REFPREFIX="$(resolve_refprefix "$REFBASE")"  # refs/hg38  OR refs/hg38.fa

  TRIM1="$WORK/${ACC}_R1.trim.fastq.gz"
  TRIM2="$WORK/${ACC}_R2.trim.fastq.gz"
  NSORT="$ALIGN/${ACC}.nsort.bam"
  FIXM="$ALIGN/${ACC}.fixmate.bam"
  SORTB="$ALIGN/${ACC}.sorted.bam"
  DEDUP="$ALIGN/${ACC}.dedup.bam"

  if [[ "$FORCE" -eq 1 ]]; then
    rm -f "$NSORT" "$FIXM" "$SORTB" "$SORTB.bai" "$DEDUP" "$DEDUP.bai" 2>/dev/null || true
  fi

  if [[ -s "$DEDUP" && -s "$DEDUP.bai" ]]; then
    echo "$ACC: dedup BAM exists -> skip alignment"
    continue
  fi

  echo "$ACC: bwa mem"
  bwa mem -t "$THREADS" "$REFPREFIX" "$TRIM1" "$TRIM2" \
  | samtools sort -n -@ "$THREADS" -o "$NSORT" -

  echo "$ACC: fixmate"
  samtools fixmate -m -@ "$THREADS" "$NSORT" "$FIXM"

  echo "$ACC: coord sort"
  samtools sort -@ "$THREADS" -o "$SORTB" "$FIXM"

  echo "$ACC: markdup"
  samtools markdup -r -@ "$THREADS" "$SORTB" "$DEDUP"
  samtools index -@ "$THREADS" "$DEDUP"

  rm -f "$NSORT" "$FIXM" "$SORTB" "$SORTB.bai" || true
done

# -------------------- 3) MACS peaks + blacklist --------------------
echo "=== [3/4] macs peaks + blacklist ==="
for ACC in $SAMPLES; do
  DEDUP="$ALIGN/${ACC}.dedup.bam"
  [[ -s "$DEDUP" ]] || { echo "Missing $DEDUP"; exit 1; }

  GSIZE="$(gsize_for "$ACC")"         # hs
  RAW="$PEAKS/${ACC}_peaks.narrowPeak"
  FILT="$PEAKS/${ACC}.peaks.filt.bed"

  if [[ "$FORCE" -eq 1 || ! -s "$RAW" ]]; then
    echo "$ACC: calling peaks"
    $MACS_BIN callpeak -t "$DEDUP" -f BAMPE -g "$GSIZE" -n "$ACC" \
      --outdir "$PEAKS" --nomodel --shift -100 --extsize 200 --keep-dup all --pvalue 1e-3 \
      >/dev/null 2>&1
  else
    echo "$ACC: peaks exist -> skip callpeak"
  fi

  # blacklist-filtered peaks
  BLK="$REFDIR/hg38.blacklist.bed"
  [[ -s "$BLK" ]] || { echo "Missing blacklist: $BLK"; exit 1; }

  if [[ "$FORCE" -eq 1 || ! -s "$FILT" ]]; then
    grep -v '^#' "$RAW" \
    | bedtools intersect -v -a - -b "$BLK" \
    > "$FILT"
    echo "$ACC: wrote $FILT"
  else
    echo "$ACC: filtered peaks exist -> skip"
  fi
done

# -------------------- 4) FRiP (fragments ∩ filtered peaks) --------------------
echo "=== [4/4] FRiP ==="
# fresh header every run (single-sample task)
echo -e "sample\tRINP\tPP_nuclear\tFRiP" > "$QC/frip.tsv"

for ACC in $SAMPLES; do
  DEDUP="$ALIGN/${ACC}.dedup.bam"
  FILT="$PEAKS/${ACC}.peaks.filt.bed"
  FRAGS="$WORK/${ACC}.frags.bed"

  # Properly paired, primary, mapped, non-dup -> BEDPE -> fragments
  samtools view -@ "$THREADS" -f 0x2 -F 0x904 -b "$DEDUP" \
  | bedtools bamtobed -bedpe -mate1 \
  | awk 'BEGIN{OFS="\t"} ($1 !~ /_/) && $2>=0 && $6>$2 {print $1,$2,$6}' \
  > "$FRAGS"

  # RINP = fragments overlapping filtered peaks
  RINP=$(bedtools intersect -u -a "$FRAGS" -b "$FILT" | wc -l | tr -d ' ')
  # PP_nuclear = count properly paired, primary, mapped, non-dup
  PP=$(samtools view -c -f 0x2 -F 0x904 "$DEDUP")
  # FRiP as fraction (0-1); write as decimal, Q&A will convert to %
  FRIP=$(python3 - <<PY
rinp=$RINP; pp=$PP
print(f"{rinp/pp:.3f}" if pp>0 else "0.000")
PY)

  printf "%s\t%s\t%s\t%s\n" "$ACC" "$RINP" "$PP" "$FRIP" >> "$QC/frip.tsv"
  echo "FRiP -> $ACC : RINP=$RINP, PP_nuclear=$PP, FRiP=$FRIP"
done

echo "Done ✅"
