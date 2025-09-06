# slsc-atac_seq analysis
## 1. Data sources
This repo demonstrates a minimal ATAC-seq workflow using one subsampled paired-end dataset (SRX26680000, FFPE Never-met SCLC tumor, aligned to hg38).

- **GEO Series:** **GSE281523** — *FFPE-ATAC and PDX-ATAC sequence data of small cell lung cancer*
- **BioProject:** **PRJNA1184329**
- **Paper:** Kawasaki K, Salehi S, Zhan YA, Chen K, *et al.* **FOXA2 promotes metastatic competence in small cell lung cancer.** *Nat Commun* 2025;16:4865.
- **Assay:** ATAC-seq (paired-end, NovaSeq 6000)

**Subset actually used for this repo:**
	•	SRX26680000 — FFPE, Never-met primary tumor (rep1)
	•	Input: one paired SRR subsampled to <1 GB each
	•	Reference: hg38

**Note:** Other SRX/SRR files from the study are listed in GEO but not used in this repo. Only the above subsampled pair is processed end-to-end in the workflow.

Subsampled FASTQs are placed in `atac_analysis/data/processed/` and used as the inputs for the pipeline. No network access is required to run the workflow.

## 2. Repository layout
```bash
sclc_atac_analysis/
├── README.md
└── atac_analysis/
    ├── data/
    │   ├── raw/          # FASTq from SRA
    │   └── processed/    # subsampled FASTQs (<1 GB each) used by the workflow
    ├── refs/             # hg38 + indexes + blacklists
    ├── align/            # BAMs
    ├── peaks/            # MACS peaks
    ├── qc/               # fastp reports + frip.tsv
    ├── results/          # (optional) extra outputs (e.g., bigWigs)
    ├── metadata.yaml
    ├── questions.yaml
    ├── answers.yaml
    └── workflow/
        ├── 01_fetch_sra.sh      # download FASTQs from SRA
        ├── 02_build_refs.sh     # build hg38 + blacklists
        ├── run_pipeline_min.sh  # trim → align → dedup → peaks → FRiP
        ├── make_qna.py          # generates questions/answers from outputs
        └── outputs/             # figures/zips, if any
```
## 3. Installation

**Requirements:**
	•	**Shell tools:** fastp, bwa (≥0.7.17), samtools, bedtools, pigz
	•	**Peaks:** macs3 or macs2 (Apple Silicon tip: python3 -m pip install macs3==3.0.3)
	•	**Python 3.9+ with PyYAML:** python3 -m pip install pyyaml
	•	fastq-dl or sra-toolkit - if you want to re-download from SRA
	•	seqtk - if you want the fetch script to auto-subsample
 
```bash
conda create -n atac-mini -y python=3.11
conda activate atac-mini
```
```bash
# tools via brew/apt/conda as you prefer; macs3 via pip works well:
python3 -m pip install macs3==3.0.3 pyyaml
```
## 4. Pre-processing / subsampling
Sub-sampled reads placed like this:
```bash
atac_analysis/data/processed/
├─ SRX26680000_R1.sub.fastq.gz
├─ SRX26680000_R2.sub.fastq.gz
```
If you want to (re)download & subsample from SRA:
```bash
cd atac_analysis/workflow
chmod +x 01_fetch_sra.sh
PROJECT="../.." ./01_fetch_sra.sh
```
	•	Prefetch (default): saves gz FASTQs to data/raw/, optional subsamples to data/processed/.
	•	Streaming option (edit script) can avoid storing .sra files locally.
## 5. Build references 
```bash
atac_analysis/refs/
  hg38.fa, hg38.{amb,ann,bwt,pac,sa}, hg38.chrom.sizes, hg38.blacklist.bed
```
To build from scratch:
```bash
cd atac_analysis/workflow
chmod +x 02_build_refs.sh
PROJECT="../.." ./02_build_refs.sh
```
## 6. Workflow
### Step 1 — Run minimal pipeline

**Purpose:** Trim/QC → Align → Deduplicate → Peak call → FRiP
**Tools:** fastp, bwa-mem, samtools, bedtools, macs3/macs2
```bash
cd atac_analysis/workflow
chmod +x run_pipeline_min.sh
PROJECT="../.." ./run_pipeline_min.sh
```
	•	Idempotent: it detects existing outputs and skips finished steps.
	•	Force recompute:
 ```bash
FORCE=1 PROJECT="../.." ./run_pipeline_min.sh
```
**Outputs:**
	•	qc/*.fastp.html, qc/*.fastp.json, qc/frip.tsv
	•	align/*.dedup.bam (+ .bai)
	•	peaks/*_peaks.narrowPeak, peaks/*.peaks.filt.bed
### Step 2 — Generate questions & answers
**Purpose:** Auto-produce questions.yaml and answers.yaml from your outputs

## 7. Primary outputs
	•	qc/frip.tsv — FRiP per sample
	•	peaks/SRX26680000_peaks.narrowPeak — raw peaks (FFPE Never-met)
	•	peaks/SRX26680000.peaks.filt.bed — blacklist-filtered peaks
	•	align/SRX*_dedup.bam — de-duplicated BAMs
	•	questions.yaml, answers.yaml — 5 Q&As from different pipeline stages
**Read count check (example):**
```bash
zcat atac_analysis/data/processed/SRX26680000_R1.sub.fastq.gz | wc -l
```

## 8. Notes:
	•	FRiP calculation uses properly paired, primary, mapped, non-dup reads only and robust fragment derivation.
	•	Peak calling uses MACS in BAMPE mode with --pvalue 1e-3, --nomodel, --shift -100, --extsize 200, and ENCODE blacklists.
