#!/usr/bin/env python3
"""
NF1 exon 31 (legacy 23a) percent spliced in (PSI) from RNA-seq junction reads.

Gives an orthogonal estimate of isoform II fractional abundance in the same
NGN2 samples measured by ddPCR, using the same formula as the GTEx analysis:

    PSI = mean(left inclusion, right inclusion) / (mean(left, right) + skipping)

Inputs (either or both, files or folders):
  * STAR junction tables  *.SJ.out.tab   (nf-core/rnaseq: star_salmon/log/)
      No extra packages needed. Recommended.
  * BAM files             *.bam (+ .bai)  (nf-core/rnaseq: star_salmon/)
      Needs pysam:  pip install pysam

Usage:
  python nf1_exon31_psi.py star_salmon/log/
  python nf1_exon31_psi.py star_salmon/*.markdup.sorted.bam --min-reads 10

Outputs (in the current folder):
  nf1_exon31_psi_samples.csv   one row per sample
  nf1_exon31_psi_by_day.csv    per line and day, and pooled across lines

Coordinates: GRCh38, NF1 on the plus strand.
Exon 31 = chr17:31,252,938-31,253,000 (63 bp).
"""
import argparse
import csv
import glob
import math
import os
import re
import sys
from collections import defaultdict

EXON_START = 31_252_938          # first base of exon 31
EXON_END = 31_253_000            # last base of exon 31
LEFT_INTRON_END = EXON_START - 1     # intron upstream of exon 31 ends here
RIGHT_INTRON_START = EXON_END + 1    # intron downstream of exon 31 starts here
REGION = (31_150_000, 31_350_000)    # window searched in BAM files
CHROMS = {"17", "chr17"}

DAY_RE = re.compile(r"[_.\-]?[Dd](?:ay|AY|Ay)?[_.\-]?(\d+)")


# ---------------------------------------------------------------- input readers
def junctions_from_sj(path):
    """STAR SJ.out.tab: col1 chrom, col2/3 intron first/last base (1-based),
    col7 uniquely mapping reads. Returns {(start, end): reads}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[0] not in CHROMS:
                continue
            s, e = int(f[1]), int(f[2])
            if REGION[0] <= s <= REGION[1]:
                out[(s, e)] = out.get((s, e), 0) + int(f[6])
    return out


def junctions_from_bam(path):
    """Count fragments (read pairs counted once) spanning each intron.
    Uniquely mapped (MAPQ 255 for STAR), primary, non-duplicate reads only."""
    try:
        import pysam
    except ImportError:
        sys.exit("Reading BAM files needs pysam: pip install pysam "
                 "(or use the SJ.out.tab files instead).")
    seen = defaultdict(set)
    with pysam.AlignmentFile(path, "rb") as bam:
        chrom = next((c for c in bam.references if c in CHROMS), None)
        if chrom is None:
            sys.exit(f"{path}: no chromosome named 17 or chr17")
        for r in bam.fetch(chrom, REGION[0] - 1, REGION[1]):
            if (r.is_unmapped or r.is_secondary or r.is_supplementary
                    or r.is_duplicate or r.mapping_quality < 255):
                continue
            pos = r.reference_start  # 0-based
            for op, length in r.cigartuples:
                if op == 3:  # N = intron
                    seen[(pos + 1, pos + length)].add(r.query_name)
                if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume reference
                    pos += length
    return {k: len(v) for k, v in seen.items()}


# ---------------------------------------------------------------- PSI
def find_coordinates(all_j):
    """Learn the exon 30 donor and exon 32 acceptor from all samples pooled,
    so that samples with zero inclusion reads still get PSI = 0."""
    tot = defaultdict(int)
    for j in all_j:
        for k, v in j.items():
            tot[k] += v
    left = {k: v for k, v in tot.items() if k[1] == LEFT_INTRON_END}
    right = {k: v for k, v in tot.items() if k[0] == RIGHT_INTRON_START}
    if left and right:
        lk, rk = max(left, key=left.get), max(right, key=right.get)
        return lk, rk, (lk[0], rk[1])
    # no inclusion reads anywhere: use the busiest junction spanning exon 31
    span = {k: v for k, v in tot.items() if k[0] < EXON_START and k[1] > EXON_END}
    if not span:
        return None
    sk = max(span, key=span.get)
    return (sk[0], LEFT_INTRON_END), (RIGHT_INTRON_START, sk[1]), sk


def exon31_counts(j, coords):
    lk, rk, sk = coords
    return {"exon30_end": lk[0] - 1, "exon32_start": rk[1] + 1,
            "left_incl": j.get(lk, 0), "right_incl": j.get(rk, 0), "skip": j.get(sk, 0)}


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def parse_name(path):
    name = os.path.basename(path)
    for suf in (".SJ.out.tab", ".markdup.sorted.bam", ".sorted.bam", ".bam"):
        if name.endswith(suf):
            name = name[: -len(suf)]
    m = DAY_RE.search(name)
    day = int(m.group(1)) if m else None
    line = name[: m.start()] if m else name
    return name, line.strip("_.-") or name, day


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="SJ.out.tab/BAM files or folders")
    ap.add_argument("--min-reads", type=int, default=10,
                    help="minimum informative reads per sample (default 10)")
    a = ap.parse_args()

    files = []
    for p in a.inputs:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, "*.SJ.out.tab")))
            files += sorted(glob.glob(os.path.join(p, "*.bam")))
        else:
            files.append(p)
    if not files:
        sys.exit("No SJ.out.tab or BAM files found.")

    loaded = []
    for f in files:
        j = junctions_from_sj(f) if f.endswith(".tab") else junctions_from_bam(f)
        loaded.append((f, j))
    coords = find_coordinates([j for _, j in loaded])
    if coords is None:
        sys.exit("No junctions spanning NF1 exon 31 found. Check chromosome naming and genome build.")

    rows = []
    for f, j in loaded:
        name, line, day = parse_name(f)
        c = exon31_counts(j, coords)
        incl = (c["left_incl"] + c["right_incl"]) / 2
        n = incl + c["skip"]
        if n == 0:
            print(f"{name}: no informative reads, skipped")
            continue
        rows.append({
            "sample": name, "line": line, "day": day,
            "left_incl": c["left_incl"], "right_incl": c["right_incl"],
            "skip": c["skip"], "informative": round(n, 1),
            "psi_pct": round(100 * incl / n, 2),
            "passes_min_reads": n >= a.min_reads,
            "exon30_end": c["exon30_end"], "exon32_start": c["exon32_start"],
        })

    if not rows:
        sys.exit("No samples had exon 31 junctions.")

    # sanity check on junction coordinates
    ends = {(r["exon30_end"], r["exon32_start"]) for r in rows}
    print(f"Exon 31: chr17:{EXON_START:,}-{EXON_END:,} ({EXON_END - EXON_START + 1} bp)")
    print("Flanking exon boundaries found (exon 30 end, exon 32 start):", ends)

    with open("nf1_exon31_psi_samples.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # per line and day, plus pooled across lines (read-weighted, Wilson 95% CI)
    groups = defaultdict(list)
    for r in rows:
        if r["passes_min_reads"]:
            groups[(r["line"], r["day"])].append(r)
            groups[("ALL", r["day"])].append(r)
    out = []
    for (line, day), g in sorted(groups.items(), key=lambda x: (x[0][1] or 0, x[0][0])):
        psis = [r["psi_pct"] for r in g]
        k = sum((r["left_incl"] + r["right_incl"]) / 2 for r in g)
        n = sum(r["informative"] for r in g)
        lo, hi = wilson(k, n)
        mean = sum(psis) / len(psis)
        sd = (math.sqrt(sum((x - mean) ** 2 for x in psis) / (len(psis) - 1))
              if len(psis) > 1 else "")
        out.append({
            "line": line, "day": day, "n_samples": len(g),
            "mean_psi_pct": round(mean, 2),
            "sd_psi_pct": round(sd, 2) if sd != "" else "",
            "pooled_psi_pct": round(100 * k / n, 2),
            "pooled_ci95_low": round(100 * lo, 2), "pooled_ci95_high": round(100 * hi, 2),
            "informative_reads": round(n, 1),
        })
    with open("nf1_exon31_psi_by_day.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    print(f"\n{'sample':28s} {'L':>6} {'R':>6} {'skip':>7} {'PSI %':>7}")
    for r in rows:
        flag = "" if r["passes_min_reads"] else "  (below min reads)"
        print(f"{r['sample']:28s} {r['left_incl']:>6} {r['right_incl']:>6} "
              f"{r['skip']:>7} {str(r['psi_pct']):>7}{flag}")
    print("\nWrote nf1_exon31_psi_samples.csv and nf1_exon31_psi_by_day.csv")


if __name__ == "__main__":
    main()
