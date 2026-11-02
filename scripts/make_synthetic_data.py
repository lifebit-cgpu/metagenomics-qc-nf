#!/usr/bin/env python3
"""Generate a small, deterministic synthetic paired-end shotgun metagenomics test set.

Everything is random sequence (no real genomes are used), so the data are tiny, licence-free
and the true answer is known exactly. Outputs (in --outdir):
  reads/<sample>_R{1,2}.fastq.gz   100 bp paired reads
  samplesheet.csv                  sample,fastq_1,fastq_2,metagenomic
  host_genome.fa                   synthetic 'host' genome (taxid 9606) for bowtie2
  kraken_library.fa                synthetic genomes with kraken:taxid headers (host + 6 bacteria + PhiX)
  taxonomy/nodes.dmp, names.dmp    minimal NCBI-style taxonomy for kraken2-build
  ground_truth.csv                 per-sample truth and expected QC flag
  ground_truth_taxa.csv            per-sample, per-taxon simulated pairs and expected share
"""
import argparse
import csv
import gzip
import os
import random

READ_LEN = 100
ERR_RATE = 0.003
SHORT_INSERT_FRAC = 0.10
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = str.maketrans("ACGTN", "TGCAN")

# taxid, parent, rank, scientific name
TAXONOMY = [
    (1, 1, "no rank", "root"),
    (131567, 1, "no rank", "cellular organisms"),
    (2, 131567, "superkingdom", "Bacteria"),
    (2759, 131567, "superkingdom", "Eukaryota"),
    (10239, 1, "superkingdom", "Viruses"),
    (816, 2, "genus", "Bacteroides"),
    (817, 816, "species", "Bacteroides fragilis"),
    (216851, 2, "genus", "Faecalibacterium"),
    (853, 216851, "species", "Faecalibacterium prausnitzii"),
    (561, 2, "genus", "Escherichia"),
    (562, 561, "species", "Escherichia coli"),
    (239934, 2, "genus", "Akkermansia"),
    (239935, 239934, "species", "Akkermansia muciniphila"),
    (1678, 2, "genus", "Bifidobacterium"),
    (1680, 1678, "species", "Bifidobacterium adolescentis"),
    (1350, 2, "genus", "Enterococcus"),
    (1351, 1350, "species", "Enterococcus faecalis"),
    (9605, 2759, "genus", "Homo"),
    (9606, 9605, "species", "Homo sapiens"),
    (1910954, 10239, "genus", "Sinsheimervirus"),
    (10847, 1910954, "species", "Escherichia virus phiX174"),
]
NAMES = {t[0]: t[3] for t in TAXONOMY}

# short code -> (taxid, synthetic genome length)
GENOMES = {
    "host":  (9606, 200000),
    "Bfrag": (817, 40000),
    "Fpra":  (853, 40000),
    "Ecoli": (562, 40000),
    "Amuc":  (239935, 40000),
    "Bado":  (1680, 40000),
    "Efae":  (1351, 40000),
    "phix":  (10847, 5386),
}

GUT = {"Bfrag": .40, "Fpra": .25, "Ecoli": .15, "Amuc": .10, "Bado": .06, "Efae": .04}
SAMPLES = [
    dict(sample="S01_typical", metagenomic=True, pairs=20000, host=.05, unknown=.05, phix=0.0, lowqual=0.0,
         taxa=GUT, expected_flag="PASS", note="typical gut sample, low host"),
    dict(sample="S02_high_host", metagenomic=True, pairs=20000, host=.60, unknown=.04, phix=0.0, lowqual=0.0,
         taxa={"Bfrag": .30, "Fpra": .30, "Ecoli": .20, "Amuc": .10, "Bado": .05, "Efae": .05},
         expected_flag="WARN", note="60% host reads (e.g. biopsy-like)"),
    dict(sample="S03_low_depth", metagenomic=True, pairs=2000, host=.10, unknown=.05, phix=0.0, lowqual=0.0,
         taxa=GUT, expected_flag="FAIL", note="below min_reads"),
    dict(sample="S04_flag_false", metagenomic=False, pairs=15000, host=.20, unknown=.05, phix=0.0, lowqual=0.0,
         taxa={"Ecoli": .50, "Efae": .30, "Bfrag": .20}, expected_flag="PASS",
         note="metagenomic=FALSE: host depletion must be skipped; host reads stay and are classified as Homo sapiens"),
    dict(sample="S05_low_quality", metagenomic=True, pairs=15000, host=.08, unknown=.05, phix=0.0, lowqual=.30,
         taxa=GUT, expected_flag="WARN", note="30% of pairs low quality, removed by fastp"),
    dict(sample="S06_phix_contam", metagenomic=True, pairs=15000, host=.05, unknown=.05, phix=.10, lowqual=0.0,
         taxa={"Fpra": .35, "Bfrag": .35, "Amuc": .20, "Bado": .10}, expected_flag="WARN",
         note="10% PhiX spike-in contamination"),
]


def revcomp(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n):
    return "".join(rng.choices("ACGT", k=n))


def allocate(total, fracs):
    raw = {k: total * v for k, v in fracs.items()}
    cnt = {k: int(raw[k]) for k in raw}
    rem = total - sum(cnt.values())
    for k in sorted(raw, key=lambda k: (raw[k] - cnt[k], k), reverse=True)[:rem]:
        cnt[k] += 1
    return cnt


def make_fragment(rng, genome, short):
    L = rng.randint(60, 95) if short else max(READ_LEN + 20, int(rng.gauss(300, 40)))
    if genome is None:
        return rand_seq(rng, L)
    L = min(L, len(genome))
    start = rng.randint(0, len(genome) - L)
    frag = genome[start:start + L]
    return revcomp(frag) if rng.random() < 0.5 else frag


def to_read(rng, seq, adapter):
    r = seq[:READ_LEN]
    if len(r) < READ_LEN:
        r = (r + adapter)[:READ_LEN]
        if len(r) < READ_LEN:
            r += rand_seq(rng, READ_LEN - len(r))
    r = list(r)
    nerr = sum(1 for _ in range(READ_LEN) if rng.random() < ERR_RATE)
    for i in rng.sample(range(READ_LEN), nerr):
        r[i] = rng.choice([b for b in "ACGT" if b != r[i]])
    return "".join(r)


def qual(rng, low):
    if low:
        return "".join(chr(33 + (rng.randint(2, 12) if rng.random() < 0.6 else rng.randint(30, 40)))
                       for _ in range(READ_LEN))
    return "".join(chr(33 + rng.randint(30, 40)) for _ in range(READ_LEN))


def write_fasta(fh, name, seq):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="synthetic_data")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    out = a.outdir
    os.makedirs(os.path.join(out, "reads"), exist_ok=True)
    os.makedirs(os.path.join(out, "taxonomy"), exist_ok=True)

    grng = random.Random(a.seed)
    genomes = {k: rand_seq(grng, L) for k, (tid, L) in GENOMES.items()}

    with open(os.path.join(out, "host_genome.fa"), "w") as fh:
        write_fasta(fh, "synthetic_host_chr1", genomes["host"])
    with open(os.path.join(out, "kraken_library.fa"), "w") as fh:
        for k, (tid, L) in GENOMES.items():
            write_fasta(fh, f"{k}|kraken:taxid|{tid}", genomes[k])
    with open(os.path.join(out, "taxonomy", "nodes.dmp"), "w") as fh:
        for tid, parent, rank, name in TAXONOMY:
            fh.write(f"{tid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")
    with open(os.path.join(out, "taxonomy", "names.dmp"), "w") as fh:
        for tid, parent, rank, name in TAXONOMY:
            fh.write(f"{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")

    truth_rows, taxa_rows, sheet_rows = [], [], []
    for idx, s in enumerate(SAMPLES):
        rng = random.Random(a.seed * 1000 + idx)
        n = s["pairs"]
        rest = 1.0 - s["host"] - s["unknown"] - s["phix"]
        fr = {"host": s["host"], "unknown": s["unknown"], "phix": s["phix"]}
        for k, w in s["taxa"].items():
            fr[k] = rest * w
        fr = {k: v for k, v in fr.items() if v > 0}
        cnt = allocate(n, fr)

        sources = [k for k, c in cnt.items() for _ in range(c)]
        rng.shuffle(sources)
        n_low = int(round(s["lowqual"] * n))
        low_idx = set(rng.sample(range(n), n_low))
        n_short = 0

        p1 = os.path.join(out, "reads", f"{s['sample']}_R1.fastq.gz")
        p2 = os.path.join(out, "reads", f"{s['sample']}_R2.fastq.gz")
        with gzip.open(p1, "wt", compresslevel=6) as f1, gzip.open(p2, "wt", compresslevel=6) as f2:
            for i, src in enumerate(sources):
                short = rng.random() < SHORT_INSERT_FRAC
                n_short += short
                frag = make_fragment(rng, None if src == "unknown" else genomes[src], short)
                r1 = to_read(rng, frag, ADAPTER_R1)
                r2 = to_read(rng, revcomp(frag), ADAPTER_R2)
                low = i in low_idx
                label = "unknown" if src == "unknown" else str(GENOMES[src][0])
                name = f"@{s['sample']}.{i + 1} src={label}"
                f1.write(f"{name}\n{r1}\n+\n{qual(rng, low)}\n")
                f2.write(f"{name}\n{r2}\n+\n{qual(rng, low)}\n")

        host = cnt.get("host", 0)
        unk = cnt.get("unknown", 0)
        phix = cnt.get("phix", 0)
        microbes = sum(c for k, c in cnt.items() if k not in ("host", "unknown"))  # bacteria + phix
        if s["metagenomic"]:
            kraken_in = n - host
            classified = microbes
        else:
            kraken_in = n
            classified = microbes + host
        truth_rows.append(dict(
            sample=s["sample"], metagenomic="TRUE" if s["metagenomic"] else "FALSE", raw_pairs=n,
            host_pairs=host, unknown_pairs=unk, phix_pairs=phix, lowqual_pairs=n_low, short_insert_pairs=n_short,
            expected_host_fraction=round(host / n, 4) if s["metagenomic"] else "NA",
            simulated_host_fraction=round(host / n, 4),
            expected_retention=round(1 - n_low / n, 4),
            expected_classified_fraction=round(classified / kraken_in, 4),
            expected_contaminant_fraction=round(phix / kraken_in, 4),
            expected_qc_flag=s["expected_flag"], note=s["note"]))
        for k, c in cnt.items():
            if k == "unknown" or (k == "host" and s["metagenomic"]):
                continue
            tid = GENOMES[k][0]
            taxa_rows.append(dict(sample=s["sample"], taxid=tid, name=NAMES[tid], simulated_pairs=c,
                                  expected_share_of_classified=round(c / classified, 4)))
        sheet_rows.append(dict(sample=s["sample"], fastq_1=f"reads/{s['sample']}_R1.fastq.gz",
                               fastq_2=f"reads/{s['sample']}_R2.fastq.gz",
                               metagenomic="TRUE" if s["metagenomic"] else "FALSE"))
        print(f"{s['sample']}: {n} pairs, host={host}, unknown={unk}, phix={phix}, lowqual={n_low}, short={n_short}")

    def dump(path, rows):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    dump(os.path.join(out, "samplesheet.csv"), sheet_rows)
    dump(os.path.join(out, "ground_truth.csv"), truth_rows)
    dump(os.path.join(out, "ground_truth_taxa.csv"), taxa_rows)
    print("done:", out)


if __name__ == "__main__":
    main()
