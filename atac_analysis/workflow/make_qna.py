#!/usr/bin/env python3
import json, subprocess, statistics, sys
from pathlib import Path
import yaml

# ---------------- paths ----------------
PROJECT = Path(__file__).resolve().parents[1]   # ../ (atac_analysis/)
QC      = PROJECT / "qc"
ALIGN   = PROJECT / "align"
PEAKS   = PROJECT / "peaks"

SAMPLE  = "SRX26680000"
FASTP_JSON = QC / f"{SAMPLE}.fastp.json"
BAM        = ALIGN / f"{SAMPLE}.dedup.bam"
NARROWPEAK = PEAKS / f"{SAMPLE}_peaks.narrowPeak"
FRIP_TSV   = QC / "frip.tsv"

# ---------------- tiny helpers ----------------
def run_cmd(cmd):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    return r.stdout.strip()

def require(p: Path, desc: str):
    if not p.exists():
        sys.exit(f"[ERROR] Missing {desc}: {p}")

# ---------------- require inputs ----------------
require(FASTP_JSON, "fastp JSON")
require(BAM, "deduplicated BAM")
require(NARROWPEAK, "MACS narrowPeak")
require(FRIP_TSV, "FRiP table")

# ---------------- Q1: % passed filter (fastp) ----------------
with open(FASTP_JSON) as fh:
    fj = json.load(fh)
p = fj["filtering_result"]["passed_filter_reads"]
f = fj["filtering_result"]["low_quality_reads"] + fj["filtering_result"]["too_many_N_reads"] + fj["filtering_result"]["too_short_reads"]
q1_percent_pass = round(p * 100.0 / (p + f), 2)

# ---------------- Q2: properly paired primary nuclear reads ----------------
# Total properly paired, primary, non-dup, non-QC-fail
pp = int(run_cmd(["samtools","view","-c","-f","0x2","-F","0x904", str(BAM)]))

# chrM properly paired (to exclude mitochondrial)
# - Include primary, non-dup, non-QC-fail like above, and select chrM only
try:
    chrM_pp = int(run_cmd(["bash","-lc", f"samtools view -f 0x2 -F 0x904 {BAM} chrM | wc -l"]))
except subprocess.CalledProcessError:
    # If no chrM contig exists in header, treat as zero
    chrM_pp = 0

q2_nuclear_pp = max(0, pp - chrM_pp)

# ---------------- Q3: chromosome with most peaks ----------------
chr_counts = {}
with open(NARROWPEAK) as f:
    for line in f:
        if not line.strip(): continue
        chrom = line.split("\t",1)[0]
        chr_counts[chrom] = chr_counts.get(chrom, 0) + 1
if not chr_counts:
    sys.exit("[ERROR] No peaks found in narrowPeak.")

q3_top_chr = max(chr_counts.items(), key=lambda kv: kv[1])[0]

# ---------------- Q4: median peak length (bp) ----------------
lengths = []
with open(NARROWPEAK) as f:
    for line in f:
        a = line.rstrip("\n").split("\t")
        if len(a) < 3: continue
        try:
            start = int(a[1]); end = int(a[2])
        except ValueError:
            continue
        if end > start >= 0:
            lengths.append(end - start)
if not lengths:
    sys.exit("[ERROR] Could not parse peak lengths from narrowPeak.")
q4_median_len = int(round(statistics.median(lengths)))

# ---------------- Q5: FRiP (%) from frip.tsv ----------------
# frip.tsv: header -> sample  reads_in_peaks  proper_pairs  frip
q5_frip_pct = None
with open(FRIP_TSV) as f:
    header = f.readline()
    for line in f:
        if not line.strip(): continue
        s, _, _, frip = line.strip().split("\t")
        if s == SAMPLE:
            q5_frip_pct = round(float(frip)*100.0, 2)
            break
if q5_frip_pct is None:
    sys.exit(f"[ERROR] SRX row not found in FRiP table: {SAMPLE}")

# ---------------- write questions.yaml ----------------
questions_yaml = {
  "task": (
    "You have been given subsampled ATAC-seq data from an FFPE primary SCLC tumor "
    f"({SAMPLE}). Using the provided paired-end FASTQ files and generated outputs, answer the following."
  ),
  "questions": [
    {
      "id": "q1",
      "stage": "quality_control",
      "text": f"After adapter removal, what percentage of reads passed quality filtering in {SAMPLE}?",
      "answer_type": "numeric_percent",
      "tolerance": 0.5,
    },
    {
      "id": "q2",
      "stage": "alignment",
      "text": f"In {SAMPLE}, how many properly paired primary reads mapped to the nuclear genome (excluding QC-fail and duplicates)?",
      "answer_type": "integer_exact",
    },
    {
      "id": "q3",
      "stage": "chromatin_accessibility",
      "text": f"Which chromosome shows the highest number of accessible regions (peak count) in {SAMPLE}?",
      "answer_type": "string_exact",
    },
    {
      "id": "q4",
      "stage": "chromatin_accessibility",
      "text": f"What is the median peak length (bp) among MACS narrowPeak calls for {SAMPLE}?",
      "answer_type": "integer_exact",
    },
    {
      "id": "q5",
      "stage": "qc",
      "text": f"What is the FRiP (Fraction of Reads in Peaks, percent) for {SAMPLE}?",
      "answer_type": "numeric_percent",
      "tolerance": 0.5,
    },
  ]
}

with open(PROJECT / "questions.yaml", "w") as out:
    yaml.safe_dump(questions_yaml, out, sort_keys=False)

# ---------------- write answers.yaml ----------------
answers_yaml = {
  "answers": {
    "q1": q1_percent_pass,          # %
    "q2": q2_nuclear_pp,            # integer
    "q3": q3_top_chr,               # string (e.g., "chr1")
    "q4": q4_median_len,            # integer bp
    "q5": q5_frip_pct,              # %
  }
}

with open(PROJECT / "answers.yaml", "w") as out:
    yaml.safe_dump(answers_yaml, out, sort_keys=False)

print(" Wrote questions.yaml and answers.yaml")
print(f"q1: {q1_percent_pass}")
print(f"q2: {q2_nuclear_pp}")
print(f"q3: {q3_top_chr}")
print(f"q4: {q4_median_len}")
print(f"q5: {q5_frip_pct}")
