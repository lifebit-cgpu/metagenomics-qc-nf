#!/usr/bin/env python3
# Build the per-sample metagenomics QC table (and, for synthetic data, a validation report against ground truth).
# Reads fastp JSON, host-depletion counts, and classifier outputs from the current directory.
import argparse, csv, json, os, re


def fnum(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def read_samples(path):
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.lower().startswith("sample,"):
                continue
            sid, meta = line.split(",")[:2]
            out.append({"sample": sid, "metagenomic": meta.strip().lower() == "true"})
    return sorted(out, key=lambda r: r["sample"])


def parse_fastp(sid):
    with open(f"{sid}.fastp.json") as fh:
        d = json.load(fh)
    b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    raw_reads = b["total_reads"]
    raw_pairs, trim_pairs = raw_reads // 2, a["total_reads"] // 2
    ad = d.get("adapter_cutting", {}).get("adapter_trimmed_reads")
    return dict(raw_pairs=raw_pairs, trimmed_pairs=trim_pairs,
                retained_frac=trim_pairs / raw_pairs if raw_pairs else None,
                q30_raw=b.get("q30_rate"), q30_trimmed=a.get("q30_rate"), gc_trimmed=a.get("gc_content"),
                mean_len_r1_trimmed=a.get("read1_mean_length"),
                duplication_rate=d.get("duplication", {}).get("rate"),
                adapter_trimmed_read_frac=(ad / raw_reads) if (ad is not None and raw_reads) else None)


def parse_host(sid):
    p = f"{sid}.host_counts.tsv"
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        r = list(csv.DictReader(fh, delimiter="\t"))[0]
    inp, non = int(r["input_pairs"]), int(r["nonhost_pairs"])
    return dict(host_pairs=inp - non, nonhost_pairs=non, host_fraction=(inp - non) / inp if inp else None)


def parse_panel(sid):
    cp, lp = f"{sid}.panel_counts.tsv", f"{sid}.panel.bowtie2.log"
    if not os.path.exists(cp):
        return None
    total = None
    with open(lp) as fh:
        for line in fh:
            m = re.match(r"\s*(\d+) reads; of these", line)
            if m:
                total = int(m.group(1)); break
    taxa = {}
    with open(cp) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = r["reference"].split("|")[0]
            taxa[t] = taxa.get(t, 0) + int(r["reads"])
    c = sum(taxa.values())
    return dict(classifier="bowtie2_panel (read 1)", reads_assessed=total, classified_reads=c,
                classified_fraction=(c / total) if total else None, taxa=taxa)


def parse_kraken(sid):
    p = f"{sid}.kraken2.report.txt"
    if not os.path.exists(p):
        return None
    unclass = root = 0
    taxa = {}
    with open(p) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            clade, rank, name = int(parts[1]), parts[-3].strip(), parts[-1].strip()
            if rank == "U":
                unclass += clade
            elif rank == "R":
                root += clade
            elif rank == "S":
                taxa[name] = clade
    total = unclass + root
    return dict(classifier="kraken2 (species)", reads_assessed=total, classified_reads=root,
                classified_fraction=(root / total) if total else None, taxa=taxa)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--min-raw-pairs", type=int, default=5000)
    ap.add_argument("--min-retained-frac", type=float, default=0.5)
    ap.add_argument("--warn-host-frac", type=float, default=0.5)
    ap.add_argument("--warn-classified-frac", type=float, default=0.5)
    ap.add_argument("--top-n", type=int, default=5)
    a = ap.parse_args()

    rows, taxa_rows, taxa_by = [], [], {}
    for s in read_samples(a.samples):
        sid = s["sample"]
        r = {"sample": sid, "metagenomic": "TRUE" if s["metagenomic"] else "FALSE"}
        r.update(parse_fastp(sid))
        h = parse_host(sid)
        r["host_depletion_run"] = h is not None
        r.update(h or dict(host_pairs=None, nonhost_pairs=None, host_fraction=None))
        c = parse_kraken(sid) or parse_panel(sid)
        if c:
            r.update({k: v for k, v in c.items() if k != "taxa"})
            ranked = sorted(c["taxa"].items(), key=lambda kv: (-kv[1], kv[0]))
            cl = c["classified_reads"] or 0
            lst = []
            for i, (t, n) in enumerate(ranked[:a.top_n], 1):
                share = n / cl if cl else None
                lst.append((t, n, share))
                taxa_rows.append(dict(sample=sid, rank=i, taxon=t, reads=n, share_of_classified=fnum(share),
                                      share_of_assessed=fnum(n / c["reads_assessed"] if c["reads_assessed"] else None)))
            taxa_by[sid] = lst
            for i in range(3):
                r[f"top{i + 1}_taxon"] = f"{lst[i][0]} ({lst[i][2]:.1%})" if i < len(lst) and lst[i][2] is not None else "NA"
        else:
            r.update(classifier="none", reads_assessed=None, classified_reads=None, classified_fraction=None,
                     top1_taxon="NA", top2_taxon="NA", top3_taxon="NA")
        fails, warns = [], []
        if r["raw_pairs"] < a.min_raw_pairs:
            fails.append(f"raw_pairs<{a.min_raw_pairs}")
        if r["retained_frac"] is not None and r["retained_frac"] < a.min_retained_frac:
            fails.append(f"retained_frac<{a.min_retained_frac}")
        if r["host_fraction"] is not None and r["host_fraction"] > a.warn_host_frac:
            warns.append(f"host_fraction>{a.warn_host_frac}")
        if r["classified_fraction"] is not None and r["classified_fraction"] < a.warn_classified_frac:
            warns.append(f"classified_fraction<{a.warn_classified_frac}")
        r["qc_flag"] = "FAIL" if fails else ("WARN" if warns else "PASS")
        r["qc_reasons"] = ";".join(fails + warns) or "none"
        rows.append(r)
        print(f"{sid}: raw={r['raw_pairs']} trimmed={r['trimmed_pairs']} host={fnum(r['host_fraction'])} "
              f"classified={fnum(r['classified_fraction'])} flag={r['qc_flag']}")

    cols = ["sample", "metagenomic", "raw_pairs", "trimmed_pairs", "retained_frac", "q30_raw", "q30_trimmed",
            "gc_trimmed", "mean_len_r1_trimmed", "duplication_rate", "adapter_trimmed_read_frac",
            "host_depletion_run", "host_pairs", "nonhost_pairs", "host_fraction", "classifier", "reads_assessed",
            "classified_reads", "classified_fraction", "top1_taxon", "top2_taxon", "top3_taxon", "qc_flag", "qc_reasons"]
    with open("qc_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([fnum(r.get(c)) for c in cols])
    with open("top_taxa.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "rank", "taxon", "reads", "share_of_classified", "share_of_assessed"])
        w.writeheader(); w.writerows(taxa_rows)
    mcols = ["metagenomic", "raw_pairs", "retained_frac", "host_fraction", "classified_fraction", "top1_taxon", "qc_flag", "qc_reasons"]
    with open("qc_summary_mqc.tsv", "w") as fh:
        fh.write("# id: 'metagenomics_qc_summary'\n# section_name: 'Metagenomics QC summary'\n")
        fh.write("# description: 'Per-sample read retention, host fraction (metagenomic=TRUE only), classified fraction, top taxon and QC flag.'\n")
        fh.write("# plot_type: 'table'\n")
        fh.write("Sample\t" + "\t".join(mcols) + "\n")
        for r in rows:
            fh.write(r["sample"] + "\t" + "\t".join(fnum(r.get(c)) for c in mcols) + "\n")

    if os.path.exists("ground_truth.csv"):
        by = {r["sample"]: r for r in rows}
        out = []
        with open("ground_truth.csv") as fh:
            for t in csv.DictReader(fh):
                s = t["sample"]; r = by.get(s)
                def add(check, exp, obs, tol, ok):
                    out.append(dict(sample=s, check=check, expected=exp, observed=obs, tolerance=tol, result="PASS" if ok else "FAIL"))
                if r is None:
                    add("sample_present", "yes", "no", "exact", False); continue
                add("raw_pairs", t["raw_pairs"], r["raw_pairs"], "exact", int(t["raw_pairs"]) == r["raw_pairs"])
                if t["metagenomic"].upper() == "TRUE":
                    e, o = float(t["true_host_fraction"]), r["host_fraction"]
                    add("host_fraction", f"{e:.4f}", fnum(o), "+/-0.02", o is not None and abs(o - e) <= 0.02)
                else:
                    o = f"{fnum(r['host_fraction'])} / {'run' if r['host_depletion_run'] else 'not run'}"
                    add("host_depletion_skipped", "NA / not run", o, "exact", r["host_fraction"] is None and not r["host_depletion_run"])
                e, o = float(t["expected_classified_fraction"]), r["classified_fraction"]
                add("classified_fraction", f"{e:.4f}", fnum(o), "+/-0.05", o is not None and abs(o - e) <= 0.05)
                shares = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in t["species_shares"].split(";")}
                top3 = t["top3"].split(";")
                obs = taxa_by.get(s, [])
                top5 = [x[0] for x in obs[:5]]
                add("top3_in_top5", ";".join(top3), ";".join(top5), "all 3 present", all(x in top5 for x in top3))
                osh = {x[0]: x[2] for x in obs}
                for sp in top3:
                    o = osh.get(sp)
                    add(f"share:{sp}", f"{shares[sp]:.4f}", fnum(o), "+/-0.05", o is not None and abs(o - shares[sp]) <= 0.05)
                add("qc_flag", t["expected_flag"], r["qc_flag"], "exact", t["expected_flag"] == r["qc_flag"])
        with open("validation_report.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["sample", "check", "expected", "observed", "tolerance", "result"])
            w.writeheader(); w.writerows(out)
        n_fail = sum(1 for x in out if x["result"] == "FAIL")
        print(f"VALIDATION: {len(out) - n_fail}/{len(out)} checks passed")


if __name__ == "__main__":
    main()
