#!/usr/bin/env python3
# Generate a small, deterministic synthetic paired-end shotgun-metagenomics dataset with a ground-truth table.
# Everything is simulated (random genomes); no real human or microbial sequence is used.
import argparse, csv, gzip, os, random

COMP = str.maketrans("ACGTN", "TGCAN")
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
READ_LEN = 100
ERROR_RATE = 0.005
HOST_LEN = 300000
SPECIES_LEN = 60000
SPECIES = [("Bacteroides_vulgatus_synth", 0.42), ("Faecalibacterium_prausnitzii_synth", 0.56),
           ("Prevotella_copri_synth", 0.45), ("Akkermansia_muciniphila_synth", 0.55),
           ("Escherichia_coli_synth", 0.51), ("Bifidobacterium_longum_synth", 0.60)]
B, F, P, A, E, L = [s[0] for s in SPECIES]
SAMPLES = [
    dict(sample="S1_typical", metagenomic=True, pairs=20000, host=0.10, unknown=0.10, adapter=0.0, lowq=0.0,
         abund={B: .35, F: .25, P: .15, A: .10, E: .10, L: .05}),
    dict(sample="S2_high_host", metagenomic=True, pairs=20000, host=0.60, unknown=0.10, adapter=0.0, lowq=0.0,
         abund={E: .40, B: .25, L: .15, F: .10, A: .05, P: .05}),
    dict(sample="S3_flag_false", metagenomic=False, pairs=20000, host=0.05, unknown=0.10, adapter=0.0, lowq=0.0,
         abund={P: .40, B: .20, F: .20, L: .10, A: .05, E: .05}),
    dict(sample="S4_low_depth", metagenomic=True, pairs=1500, host=0.10, unknown=0.10, adapter=0.0, lowq=0.0,
         abund={A: .35, B: .25, F: .20, P: .10, E: .05, L: .05}),
    dict(sample="S5_adapters_lowq", metagenomic=True, pairs=20000, host=0.20, unknown=0.10, adapter=0.30, lowq=0.30,
         abund={F: .30, B: .25, L: .20, P: .10, A: .10, E: .05}),
]


def rc(s):
    return s.translate(COMP)[::-1]


def make_genome(rng, n, gc):
    at = (1 - gc) / 2
    return "".join(rng.choices("ACGT", weights=[at, gc / 2, gc / 2, at], k=n))


def add_errors(rng, seq, rate):
    s = list(seq)
    i = -1
    while True:
        i += 1 + int(rng.expovariate(rate))
        if i >= len(s):
            break
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def allocate(n, weights):
    # largest-remainder allocation of n items to weights (dict name -> weight)
    raw = {k: n * w for k, w in weights.items()}
    out = {k: int(v) for k, v in raw.items()}
    rest = n - sum(out.values())
    for k in sorted(raw, key=lambda k: (raw[k] - out[k], k), reverse=True)[:rest]:
        out[k] += 1
    return out


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min-raw-pairs", type=int, default=5000)
    ap.add_argument("--warn-host-frac", type=float, default=0.5)
    ap.add_argument("--warn-classified-frac", type=float, default=0.5)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    os.makedirs(os.path.join(a.outdir, "reads"), exist_ok=True)

    host = make_genome(rng, HOST_LEN, 0.41)
    genomes = {name: make_genome(rng, SPECIES_LEN, gc) for name, gc in SPECIES}
    write_fasta(os.path.join(a.outdir, "host_synthetic.fa"), [("host_synthetic_chr1", host)])
    write_fasta(os.path.join(a.outdir, "panel_synthetic.fa"), [(f"{n}|contig1", genomes[n]) for n, _ in SPECIES])

    sheet, truth = [], []
    for s in SAMPLES:
        n = s["pairs"]
        nh = round(n * s["host"])
        nu = round(n * s["unknown"])
        nb = n - nh - nu
        sp_counts = allocate(nb, s["abund"])
        sources = ["host"] * nh + ["unknown"] * nu
        for sp, c in sp_counts.items():
            sources += [sp] * c
        rng.shuffle(sources)
        f1 = f"{s['sample']}_R1.fastq.gz"
        f2 = f"{s['sample']}_R2.fastq.gz"
        with gzip.open(os.path.join(a.outdir, "reads", f1), "wt", compresslevel=6) as o1, \
             gzip.open(os.path.join(a.outdir, "reads", f2), "wt", compresslevel=6) as o2:
            for i, src in enumerate(sources):
                short = rng.random() < s["adapter"]
                ins = rng.randint(60, 95) if short else max(150, min(600, int(rng.gauss(300, 40))))
                if src == "unknown":
                    frag = make_genome(rng, ins, 0.50)
                else:
                    g = host if src == "host" else genomes[src]
                    st = rng.randrange(0, len(g) - ins)
                    frag = g[st:st + ins]
                    if rng.random() < 0.5:
                        frag = rc(frag)
                if len(frag) >= READ_LEN:
                    r1, r2 = frag[:READ_LEN], rc(frag)[:READ_LEN]
                else:
                    pad = "".join(rng.choices("ACGT", k=READ_LEN))
                    r1 = (frag + ADAPTER_R1 + pad)[:READ_LEN]
                    r2 = (rc(frag) + ADAPTER_R2 + pad)[:READ_LEN]
                r1, r2 = add_errors(rng, r1, ERROR_RATE), add_errors(rng, r2, ERROR_RATE)
                q1 = q2 = "I" * READ_LEN
                if rng.random() < s["lowq"]:
                    tail = 20
                    r1 = r1[:-tail] + "".join(rng.choices("ACGT", k=tail))
                    r2 = r2[:-tail] + "".join(rng.choices("ACGT", k=tail))
                    q1 = q2 = "I" * (READ_LEN - tail) + "#" * tail
                name = f"{s['sample']}_{i + 1:06d}_src-{src}"
                o1.write(f"@{name} 1:N:0:1\n{r1}\n+\n{q1}\n")
                o2.write(f"@{name} 2:N:0:1\n{r2}\n+\n{q2}\n")
        meta = "TRUE" if s["metagenomic"] else "FALSE"
        sheet.append(dict(sample=s["sample"], fastq_1=f1, fastq_2=f2, metagenomic=meta))
        host_frac = nh / n
        exp_class = nb / (n - nh) if s["metagenomic"] else nb / n
        shares = {sp: c / nb for sp, c in sp_counts.items()}
        top3 = sorted(shares, key=lambda k: (-shares[k], k))[:3]
        if n < a.min_raw_pairs:
            flag = "FAIL"
        elif (s["metagenomic"] and host_frac > a.warn_host_frac) or exp_class < a.warn_classified_frac:
            flag = "WARN"
        else:
            flag = "PASS"
        truth.append(dict(sample=s["sample"], metagenomic=meta, raw_pairs=n, host_pairs=nh, unknown_pairs=nu,
                          bacterial_pairs=nb, true_host_fraction=f"{host_frac:.4f}",
                          expected_classified_fraction=f"{exp_class:.4f}",
                          species_shares=";".join(f"{k}={v:.4f}" for k, v in sorted(shares.items())),
                          top3=";".join(top3), adapter_pair_frac=s["adapter"], lowq_pair_frac=s["lowq"],
                          expected_flag=flag))
        print(f"{s['sample']}: {n} pairs, host {nh}, unknown {nu}, bacterial {nb}, expected flag {flag}")

    with open(os.path.join(a.outdir, "samplesheet.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2", "metagenomic"])
        w.writeheader(); w.writerows(sheet)
    with open(os.path.join(a.outdir, "ground_truth.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(truth[0].keys()))
        w.writeheader(); w.writerows(truth)


if __name__ == "__main__":
    main()
