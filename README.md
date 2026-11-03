# metagenomics-qc-nf

Nextflow (DSL2) QC pipeline for vendor-delivered paired-end shotgun metagenomics FASTQ (for example from CosmosID, whose deliverables include taxonomic profiles but no QC metrics).

## What it does

| Step | Tool | Output |
|---|---|---|
| Raw read QC | FastQC 0.12.1 | `fastqc_raw/` |
| Adapter and quality trimming | fastp 0.23.4 (`--detect_adapter_for_pe --cut_tail`, min length 50) | `fastp/` |
| Host-read detection and depletion, **only for samples with `metagenomic = TRUE`** | Bowtie2 2.5.4 (`--sensitive`, `--un-conc-gz`) | `host_depletion/` (logs, per-sample host counts) |
| Classified fraction and top taxa | Kraken2 2.1.3 (`--kraken2_db`) **or** Bowtie2 against a reference panel (`--classifier_fasta`, headers `>taxon|contig`) | `classification/` |
| Per-sample QC table with PASS/WARN/FAIL | `bin/qc_summary.py` | `qc_summary/qc_summary.csv`, `top_taxa.csv` |
| Aggregate report | MultiQC 1.25.1 (FastQC, fastp, Bowtie2, custom QC table) | `multiqc/multiqc_report.html` |

Per-sample metrics: raw and trimmed read pairs, retained fraction, Q30 before/after, GC, duplication rate, adapter-trimmed read fraction, host pairs and **host fraction** (NA when `metagenomic = FALSE`), **classified fraction**, **top taxa**, QC flag with reasons.

QC flag rules (all thresholds are parameters): FAIL if raw pairs < `min_raw_pairs` (5000) or retained fraction < `min_retained_frac` (0.5); WARN if host fraction > `warn_host_frac` (0.5) or classified fraction < `warn_classified_frac` (0.5); otherwise PASS.

## Input

`--input samplesheet.csv`:

```
sample,fastq_1,fastq_2,metagenomic
STOOL_001,s3://bucket/STOOL_001_R1.fastq.gz,s3://bucket/STOOL_001_R2.fastq.gz,TRUE
```

`metagenomic` = TRUE/FALSE controls whether host detection and depletion run for that sample.

References: `--host_index <dir with a bowtie2 index>` (recommended for human, e.g. a pre-built GRCh38/CHM13 index) or `--host_fasta <fasta>` (index is built at run time). Classification: `--kraken2_db <dir>` or `--classifier_fasta <fasta>`; if neither is given, classified fraction and top taxa are NA.

## Test with synthetic data

```
nextflow run lifebit-cgpu/metagenomics-qc-nf -r batch-v1 -profile test,docker
```

`-profile test` generates a deterministic synthetic dataset at run time (`bin/make_test_data.py`, seed 42): a random 300 kb "host" genome, six random 60 kb "species" genomes, and five paired-end samples (2x100 bp) designed to exercise every branch:

| Sample | Pairs | Host | metagenomic | Designed to test | Expected flag |
|---|---|---|---|---|---|
| S1_typical | 20,000 | 10% | TRUE | normal sample | PASS |
| S2_high_host | 20,000 | 60% | TRUE | high host fraction | WARN |
| S3_flag_false | 20,000 | 5% | FALSE | host depletion skipped | PASS |
| S4_low_depth | 1,500 | 10% | TRUE | low depth | FAIL |
| S5_adapters_lowq | 20,000 | 20% | TRUE | adapter read-through (30% of pairs) and low-quality tails (30%) | PASS |

Each sample also contains 10% random "unknown" reads that match no reference. Read names carry the true source. `ground_truth.csv` records the true counts, and `QC_SUMMARY` writes `qc_summary/validation_report.csv`, which compares pipeline output with the truth using fixed tolerances: raw pairs exact; host fraction within 0.02; classified fraction within 0.05; true top-3 species within the top 5; each top-3 share within 0.05; QC flag exact.

All sequence is simulated; no real human or microbial data are included.

## Profiles

- `test`: synthetic data (above)
- `docker`: run tasks in the pinned biocontainers
- `local`: no containers, small resources (tools must be on PATH)

On CloudOS / AWS Batch the per-process containers in `nextflow.config` are used automatically.

## Limitations

- The bowtie2-panel classifier in the test profile counts read 1 only and only knows the panel genomes; for real data use a Kraken2 database (`--kraken2_db`). The Kraken2 path has only been exercised in stub mode.
- Host fraction counts a pair as host only if it aligns concordantly to the host reference.
