# metagenomics-qc-nf

QC pipeline (Nextflow DSL2) for vendor-delivered paired-end shotgun metagenomics FASTQ, e.g. CosmosID deliveries,
where the vendor supplies taxonomic / gut-brain / gut-metabolic module results and raw FASTQ but **no QC metrics**.

**Status: v0.1.0, not yet validated.** It has not completed a run yet; see "Validation status" below.

## What it does

```
raw FASTQ ─► FastQC (raw)
         └─► fastp (adapter + quality trimming) ─┬─ metagenomic=TRUE  ─► bowtie2 host depletion ─┐
                                                 └─ metagenomic=FALSE ─  (skipped) ───────────────┤
                                                                         Kraken2 (optional) ◄─────┘
                                           QC summary (per sample) + MultiQC report
```

Per-sample outputs (`results/qc_summary/qc_summary.csv`):

| column | meaning |
|---|---|
| raw_pairs, trimmed_pairs, retention | read pairs before/after fastp and the fraction kept |
| adapter_trimmed_reads, q30_rate_raw, gc_content_raw | from fastp |
| host_depletion | `performed`, `skipped (metagenomic=FALSE)` or `skipped (no host reference)` |
| host_pairs, host_fraction | pairs aligning to the host (bowtie2) / trimmed pairs. `NA` when depletion is skipped |
| classified_fraction | Kraken2 classified pairs / pairs given to Kraken2 (non-host pairs when depleted) |
| contaminant_fraction | pairs assigned to `--contaminant_taxids` (default PhiX 10847) / Kraken2 input |
| top_taxa | top `--top_n` species with % of classified pairs (full list in `top_taxa.csv`) |
| qc_flag, qc_reasons | FAIL if raw_pairs < `min_reads`; WARN if retention < `min_retention`, host_fraction > `max_host_fraction` or contaminant_fraction > `max_contaminant_fraction` |

MultiQC (`results/multiqc/multiqc_report.html`) combines FastQC, fastp, bowtie2, Kraken2 and two custom sections:
the QC summary table and a read-fate bar chart.

## Running

### Built-in synthetic test data
```
nextflow run lifebit-cgpu/metagenomics-qc-nf -r main -profile test,docker
```
On CloudOS: import the repository as a Nextflow pipeline and launch with profile `test`. CloudOS supplies the
AWS Batch executor and Docker settings.

`-profile test` runs `scripts/make_synthetic_data.py` (seed 42) inside the pipeline. It builds random-sequence genomes
(a 200 kb synthetic host with taxid 9606, six 40 kb gut bacteria and a 5.4 kb PhiX), a minimal NCBI-style taxonomy and a
Kraken2 library. It then simulates six samples of 100 bp read pairs (0.3% substitution errors, 10% short inserts with adapter read-through):

| sample | metagenomic | pairs | design | expected QC flag |
|---|---|---|---|---|
| S01_typical | TRUE | 20,000 | 5% host | PASS |
| S02_high_host | TRUE | 20,000 | 60% host | WARN (host) |
| S03_low_depth | TRUE | 2,000 | 10% host | FAIL (depth) |
| S04_flag_false | FALSE | 15,000 | 20% host, depletion must be skipped | PASS |
| S05_low_quality | TRUE | 15,000 | 30% low-quality pairs | WARN (retention) |
| S06_phix_contam | TRUE | 15,000 | 10% PhiX | WARN (contamination) |

The exact truth is published to `results/test_data/ground_truth.csv` and `ground_truth_taxa.csv`.
The generator can also be run on its own: `python3 scripts/make_synthetic_data.py --outdir synthetic_data --seed 42`.

### Real data
```
nextflow run lifebit-cgpu/metagenomics-qc-nf -r main -profile docker \
  --input samplesheet.csv \
  --host_fasta s3://.../GRCh38.fa        # or --bowtie2_index s3://.../bt2_dir/
  --kraken2_db s3://.../k2_standard/     # optional; omit for QC without classification
  --outdir results
```
`samplesheet.csv`:
```
sample,fastq_1,fastq_2,metagenomic
P001,s3://bucket/P001_R1.fastq.gz,s3://bucket/P001_R2.fastq.gz,TRUE
P002,s3://bucket/P002_R1.fastq.gz,s3://bucket/P002_R2.fastq.gz,FALSE
```
Relative paths are resolved against the samplesheet's folder. Only paired-end input is supported.
For a full human index or Kraken2 standard database, raise memory (for example `--kraken2_memory '80 GB'`).

## Parameters
See `nextflow.config`. Main ones: `fastp_qualified_quality` 15, `fastp_unqualified_percent` 40,
`fastp_min_length` 50, `min_reads` 5000, `max_host_fraction` 0.5, `min_retention` 0.8,
`max_contaminant_fraction` 0.01, `contaminant_taxids` '10847', `top_n` 5.

## Software (containers)
FastQC 0.12.1, fastp 0.23.4, Bowtie2 2.5.4, Kraken2 2.1.3, MultiQC 1.25.1 (quay.io/biocontainers; tags checked
on quay.io on 2026-11-02), python:3.12-slim for the generator and summary scripts. Nextflow >= 22.10 (DSL2).

## Notes and limitations
- Host fraction is "pairs not left unaligned by bowtie2 `--un-conc`" over trimmed pairs. Pairs where only one mate
  aligns count as host.
- Classified fraction depends entirely on the Kraken2 database used. CosmosID's own classifier and database differ,
  so these numbers are a QC indication, not a replacement for the vendor profile.
- The synthetic data are random sequence. They test the pipeline logic and the metric calculations, not
  real-world classification accuracy.

## Validation status
Not yet run end to end. Container tags were checked on quay.io, but no pipeline execution has completed.
