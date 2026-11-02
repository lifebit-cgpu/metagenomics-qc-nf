#!/usr/bin/env python3
"""Combine fastp, bowtie2 host-depletion and Kraken2 outputs into a per-sample QC table.

Outputs: qc_summary.csv, top_taxa.csv, qc_summary_mqc.tsv and read_fate_mqc.tsv (MultiQC custom content).
"""
import argparse
import csv
import glob
import json
import os


def fmt(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def read_kraken(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rows.append(dict(clade=int(f[1]), rank=f[3].strip(), taxid=f[4].strip(), name=f[5].strip()))
    return rows


def find_one(d, sample, suffix):
    hits = sorted(glob.glob(os.path.join(d, f"{sample}{suffix}")))
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--fastp_dir", default="fastp")
    ap.add_argument("--host_dir", default="host")
    ap.add_argument("--kraken_dir", default="kraken")
    ap.add_argument("--min_reads", type=int, default=5000)
    ap.add_argument("--max_host_fraction", type=float, default=0.5)
    ap.add_argument("--min_retention", type=float, default=0.8)
    ap.add_argument("--max_contaminant_fraction", type=float, default=0.01)
    ap.add_argument("--contaminant_taxids", default="10847")
    ap.add_argument("--top_n", type=int, default=5)
    a = ap.parse_args()
    contam = {t.strip() for t in str(a.contaminant_taxids).split(",") if t.strip()}

    samples = []
    with open(a.samples) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("sample,"):
                continue
            s, mg = line.rsplit(",", 1)
            samples.append((s, mg.strip().lower() == "true"))
    samples.sort()

    rows, taxa_rows, fate_rows = [], [], []
    for s, mg in samples:
        r = dict(sample=s, metagenomic="TRUE" if mg else "FALSE")
        fj = find_one(a.fastp_dir, s, ".fastp.json")
        raw = trimmed = None
        if fj:
            with open(fj) as fh:
                js = json.load(fh)
            raw = js["summary"]["before_filtering"]["total_reads"] // 2
            trimmed = js["summary"]["after_filtering"]["total_reads"] // 2
            r["q30_rate_raw"] = js["summary"]["before_filtering"].get("q30_rate")
            r["gc_content_raw"] = js["summary"]["before_filtering"].get("gc_content")
            r["adapter_trimmed_reads"] = js.get("adapter_cutting", {}).get("adapter_trimmed_reads")
        r["raw_pairs"] = raw
        r["trimmed_pairs"] = trimmed
        r["retention"] = (trimmed / raw) if raw else None

        hs = find_one(a.host_dir, s, ".host_stats.tsv")
        host_pairs = nonhost = None
        if mg and hs:
            with open(hs) as fh:
                rec = list(csv.DictReader(fh, delimiter="\t"))[0]
            inp, nonhost = int(rec["input_pairs"]), int(rec["nonhost_pairs"])
            host_pairs = inp - nonhost
            r["host_depletion"] = "performed"
        else:
            r["host_depletion"] = "skipped (metagenomic=FALSE)" if not mg else "skipped (no host reference)"
        r["host_pairs"] = host_pairs
        r["host_fraction"] = (host_pairs / trimmed) if (host_pairs is not None and trimmed) else None
        r["nonhost_pairs"] = nonhost

        kr = find_one(a.kraken_dir, s, ".kraken2.report.txt")
        r["kraken_input_pairs"] = r["classified_fraction"] = r["contaminant_fraction"] = None
        top = []
        classified = unclassified = None
        if kr:
            rep = read_kraken(kr)
            unclassified = sum(x["clade"] for x in rep if x["rank"] == "U")
            classified = sum(x["clade"] for x in rep if x["taxid"] == "1")
            total = classified + unclassified
            r["kraken_input_pairs"] = total
            r["classified_fraction"] = classified / total if total else None
            r["contaminant_fraction"] = (sum(x["clade"] for x in rep if x["taxid"] in contam) / total) if total else None
            sp = sorted([x for x in rep if x["rank"] == "S" and x["clade"] > 0], key=lambda x: -x["clade"])
            for rank_i, x in enumerate(sp, 1):
                share = x["clade"] / classified if classified else 0.0
                taxa_rows.append(dict(sample=s, rank=rank_i, taxid=x["taxid"], name=x["name"], pairs=x["clade"],
                                      share_of_classified=round(share, 4),
                                      is_contaminant=x["taxid"] in contam))
                if rank_i <= a.top_n:
                    top.append(f"{x['name']} ({100 * share:.1f}%)")
        r["top_taxa"] = "; ".join(top) if top else "NA"

        reasons, flag = [], "PASS"
        if raw is None or raw < a.min_reads:
            flag = "FAIL"
            reasons.append(f"raw_pairs<{a.min_reads}")
        if r["retention"] is not None and r["retention"] < a.min_retention:
            reasons.append(f"retention<{a.min_retention}")
        if r["host_fraction"] is not None and r["host_fraction"] > a.max_host_fraction:
            reasons.append(f"host_fraction>{a.max_host_fraction}")
        if r["contaminant_fraction"] is not None and r["contaminant_fraction"] > a.max_contaminant_fraction:
            reasons.append(f"contaminant_fraction>{a.max_contaminant_fraction}")
        if flag != "FAIL" and reasons:
            flag = "WARN"
        r["qc_flag"] = flag
        r["qc_reasons"] = ";".join(reasons) if reasons else "none"
        rows.append(r)

        if raw is not None and trimmed is not None:
            fate = dict(sample=s, removed_by_trimming=raw - trimmed, host=host_pairs or 0)
            if classified is not None:
                fate["classified"] = classified
                fate["unclassified"] = unclassified
            else:
                fate["not_classified"] = trimmed - (host_pairs or 0)
            fate_rows.append(fate)

    cols = ["sample", "metagenomic", "raw_pairs", "trimmed_pairs", "retention", "adapter_trimmed_reads",
            "q30_rate_raw", "gc_content_raw", "host_depletion", "host_pairs", "host_fraction", "nonhost_pairs",
            "kraken_input_pairs", "classified_fraction", "contaminant_fraction", "top_taxa", "qc_flag", "qc_reasons"]
    with open("qc_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([fmt(r.get(c)) for c in cols])

    tcols = ["sample", "rank", "taxid", "name", "pairs", "share_of_classified", "is_contaminant"]
    with open("top_taxa.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=tcols)
        w.writeheader()
        w.writerows(taxa_rows)

    mcols = ["metagenomic", "raw_pairs", "retention", "host_depletion", "host_fraction", "classified_fraction",
             "contaminant_fraction", "top_taxa", "qc_flag", "qc_reasons"]
    with open("qc_summary_mqc.tsv", "w") as fh:
        fh.write("# id: 'qc_summary'\n# section_name: 'Lifebit QC summary'\n"
                 "# description: 'Per-sample QC built from raw FASTQ: read counts, trimming retention, host fraction "
                 "(of trimmed pairs; only when metagenomic=TRUE), Kraken2 classified fraction, contaminant fraction and top species.'\n"
                 "# plot_type: 'table'\n")
        fh.write("Sample\t" + "\t".join(mcols) + "\n")
        for r in rows:
            fh.write(r["sample"] + "\t" + "\t".join(fmt(r.get(c)) for c in mcols) + "\n")

    fcats = ["removed_by_trimming", "host", "classified", "unclassified", "not_classified"]
    with open("read_fate_mqc.tsv", "w") as fh:
        fh.write("# id: 'read_fate'\n# section_name: 'Read fate'\n"
                 "# description: 'Read pairs removed by trimming, identified as host, and classified/unclassified by Kraken2.'\n"
                 "# plot_type: 'bargraph'\n# pconfig:\n#     id: 'read_fate_plot'\n#     title: 'Read fate per sample'\n#     ylab: 'Read pairs'\n")
        fh.write("Sample\t" + "\t".join(fcats) + "\n")
        for f in fate_rows:
            fh.write(f["sample"] + "\t" + "\t".join(str(f.get(c, 0)) for c in fcats) + "\n")

    for r in rows:
        print(r["sample"], r["qc_flag"], "raw", r["raw_pairs"], "host_frac", fmt(r["host_fraction"]),
              "classified", fmt(r["classified_fraction"]), "| ", r["top_taxa"])


if __name__ == "__main__":
    main()
