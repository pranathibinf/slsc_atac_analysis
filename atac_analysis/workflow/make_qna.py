#!/usr/bin/env python3
import json, subprocess, sys
from pathlib import Path
import yaml

PROJECT = Path(__file__).resolve().parents[1]  # ../ (atac_analysis/)
QC_DIR   = PROJECT / "qc"
ALIGN_DIR= PROJECT / "align"
PEAKS_DIR= PROJECT / "peaks"

SAMPLE_FASTP_FOR_Q1 = "SRX26680004"
SAMPLE_ALIGN_FOR_Q2 = "SRX26680000"
SAMPLE_PEAKS_FOR_Q4 = "SRX26680000"
SAMPLE_FRIP_FOR_Q5  = "SRX26680000"

def read_fastp_json(sample: str):
    jf = QC_DIR / f"{sample}.fastp.json"
    if not jf.exists(): sys.exit(f"[ERROR] fastp JSON missing: {jf}")
    return json.loads(jf.read_text())

def percent_passed_from_fastp(fp):
    if "filtering_result" in fp:
        fr = fp["filtering_result"]
        passed = fr["passed_filter_reads"]
        failed = fr["low_quality_reads"] + fr["too_many_N_reads"] + fr["too_short_reads"]
        total = passed + failed
        return round(passed * 100.0 / total, 2) if total else 0.0
    before = fp["summary"]["before_filtering"]["total_reads"]
    after  = fp["summary"]["after_filtering"]["total_reads"]
    return round(after * 100.0 / before, 2) if before else 0.0

def duplication_percent_from_fastp(fp):
    frac = (fp.get("duplication",{}) or {}).get("rate", fp.get("summary",{}).get("duplication"))
    if frac is None: return None
    return round(float(frac) * 100.0, 2)

def count_proper_pairs(bam: Path):
    if not bam.exists(): sys.exit(f"[ERROR] BAM missing: {bam}")
    out = subprocess.check_output(["samtools","view","-c","-f","0x2","-F","0x904",str(bam)], text=True).strip()
    return int(out)

def count_peaks(sample: str):
    np = PEAKS_DIR / f"{sample}_peaks.narrowPeak"
    if not np.exists():
        nps = sorted(PEAKS_DIR.glob("*.narrowPeak"))
        if not nps: sys.exit("[ERROR] No .narrowPeak files in peaks/")
        np = nps[0]
    return sum(1 for _ in np.open())

def read_frip(sample: str):
    tsv = QC_DIR / "frip.tsv"
    if not tsv.exists(): sys.exit(f"[ERROR] FRiP table missing: {tsv}")
    with tsv.open() as f:
        f.readline()
        for line in f:
            fields = line.strip().split()
            if fields and fields[0] == sample:
                return round(float(fields[-1]) * 100.0, 2)
    sys.exit(f"[ERROR] Sample {sample} not in {tsv}")

questions = {
  "task": (
    "You have been given subsampled ATAC-seq data from SCLC (FFPE Never-met, "
    "FFPE Met-associated, and FOXA2+ PDX). Using the provided paired-end FASTQ files and "
    "generated outputs, answer the following."
  ),
  "questions": [
    {"id":"q1","stage":"quality_control",
     "text":f"After adapter removal, what percentage of reads passed filtering in sample {SAMPLE_FASTP_FOR_Q1}?",
     "answer_type":"numeric_percent","tolerance":0.5},
    {"id":"q2","stage":"alignment",
     "text":f"How many properly paired primary reads are mapped (excluding QC-fail and duplicates) in sample {SAMPLE_ALIGN_FOR_Q2}?",
     "answer_type":"integer_exact"},
    {"id":"q3","stage":"quality_control",
     "text":f"What is the duplication rate (percent) reported by fastp for sample {SAMPLE_FASTP_FOR_Q1}?",
     "answer_type":"numeric_percent","tolerance":0.1},
    {"id":"q4","stage":"peak_calling",
     "text":f"How many significant peaks were called by MACS (narrowPeak) for sample {SAMPLE_PEAKS_FOR_Q4}?",
     "answer_type":"integer_exact"},
    {"id":"q5","stage":"qc",
     "text":f"What is the FRiP (Fraction of Reads in Peaks, percent) for sample {SAMPLE_FRIP_FOR_Q5}?",
     "answer_type":"numeric_percent","tolerance":0.5}
  ]
}

fp = read_fastp_json(SAMPLE_FASTP_FOR_Q1)
answers = {
  "q1": percent_passed_from_fastp(fp),
  "q3": duplication_percent_from_fastp(fp),
  "q2": count_proper_pairs(ALIGN_DIR / f"{SAMPLE_ALIGN_FOR_Q2}.dedup.bam"),
  "q4": count_peaks(SAMPLE_PEAKS_FOR_Q4),
  "q5": read_frip(SAMPLE_FRIP_FOR_Q5)
}
if answers["q3"] is None: sys.exit("[ERROR] Could not parse duplication rate from fastp JSON.")

Path(PROJECT, "questions.yaml").write_text(yaml.safe_dump(questions, sort_keys=False, width=80))
Path(PROJECT, "answers.yaml").write_text(yaml.safe_dump({"answers": answers}, sort_keys=False, width=80))

print("✅ Wrote questions.yaml and answers.yaml")
for k, v in answers.items():
    print(f"{k}: {v}")

