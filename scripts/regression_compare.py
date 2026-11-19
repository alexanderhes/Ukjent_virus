#!/usr/bin/env python3
"""Compare a pipeline results folder with a baseline results folder.

Compared:
  overview/all_samples_overview.tsv   rows keyed by (sample_ID, virus_name)
  validation/genome_summary.tsv       rows keyed by (sample, group_id)
  validation/*/genomes/*.fasta        sequences, per file name

Text columns must match exactly. Numeric columns may differ by --rel-tol
(relative) or --abs-tol (absolute), whichever is larger; this absorbs the
known +-1 read pair from FASTP_DEDUP between runs. Genomes are reported as
identical, or with the number of differing positions / length change.

Output never leaves the machine, but it contains sample names: do not paste
it into commits or anything online.

Usage: regression_compare.py <baseline_dir> <new_dir> [--rel-tol 0.01] [--abs-tol 2]
Exit status 1 if anything differs beyond tolerance (genome differences included).
"""
import argparse
import csv
import glob
import os
import sys

TABLES = [
    ("overview/all_samples_overview.tsv", ("sample_ID", "virus_name")),
    ("validation/genome_summary.tsv", ("sample", "group_id")),
]


def read_table(path, key):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {tuple(r[k] for k in key): r for r in rows}, (rows[0].keys() if rows else [])


def as_num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def compare_tables(base, new, rel, abs_, out):
    fails = 0
    for rel_path, key in TABLES:
        bp, np_ = os.path.join(base, rel_path), os.path.join(new, rel_path)
        if not os.path.exists(bp) and not os.path.exists(np_):
            continue
        if not os.path.exists(bp) or not os.path.exists(np_):
            out.append(f"FAIL {rel_path}: present in only one folder")
            fails += 1
            continue
        b, bcols = read_table(bp, key)
        n, ncols = read_table(np_, key)
        if list(bcols) != list(ncols):
            out.append(f"FAIL {rel_path}: columns differ "
                       f"(only baseline: {sorted(set(bcols) - set(ncols))}, "
                       f"only new: {sorted(set(ncols) - set(bcols))})")
            fails += 1
        for k in sorted(set(b) - set(n)):
            out.append(f"FAIL {rel_path}: row missing in new run: {k}")
            fails += 1
        for k in sorted(set(n) - set(b)):
            out.append(f"FAIL {rel_path}: new row: {k}")
            fails += 1
        for k in sorted(set(b) & set(n)):
            for col in bcols:
                if col not in n[k]:
                    continue
                bv, nv = b[k][col], n[k][col]
                if bv == nv:
                    continue
                bn, nn = as_num(bv), as_num(nv)
                if bn is not None and nn is not None:
                    tol = max(abs_, rel * abs(bn))
                    if abs(bn - nn) <= tol:
                        out.append(f"ok   {rel_path}: {k} {col}: {bv} -> {nv} (within tolerance)")
                        continue
                out.append(f"FAIL {rel_path}: {k} {col}: {bv} -> {nv}")
                fails += 1
    return fails


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def compare_genomes(base, new, out):
    fails = 0
    pattern = os.path.join("validation", "*", "genomes", "*.fasta")
    bf = {os.path.relpath(p, base) for p in glob.glob(os.path.join(base, pattern))}
    nf = {os.path.relpath(p, new) for p in glob.glob(os.path.join(new, pattern))}
    for p in sorted(bf - nf):
        out.append(f"FAIL genome file missing in new run: {p}")
        fails += 1
    for p in sorted(nf - bf):
        out.append(f"FAIL new genome file: {p}")
        fails += 1
    for p in sorted(bf & nf):
        b, n = read_fasta(os.path.join(base, p)), read_fasta(os.path.join(new, p))
        if b == n:
            continue
        for name in sorted(set(b) | set(n)):
            bs, ns = b.get(name), n.get(name)
            if bs == ns:
                continue
            if bs is None or ns is None:
                out.append(f"FAIL {p}: record {name} only in {'new' if bs is None else 'baseline'}")
            elif len(bs) == len(ns):
                d = sum(1 for x, y in zip(bs, ns) if x != y)
                out.append(f"FAIL {p}: {name}: {d} of {len(bs)} positions differ")
            else:
                out.append(f"FAIL {p}: {name}: length {len(bs)} -> {len(ns)}")
            fails += 1
    return fails


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("baseline")
    ap.add_argument("new")
    ap.add_argument("--rel-tol", type=float, default=0.01)
    ap.add_argument("--abs-tol", type=float, default=2)
    a = ap.parse_args()
    out = []
    fails = compare_tables(a.baseline, a.new, a.rel_tol, a.abs_tol, out)
    fails += compare_genomes(a.baseline, a.new, out)
    for line in out:
        print(line)
    print(f"\nregression: {'PASS' if fails == 0 else f'FAIL ({fails} differences)'}"
          f" (baseline {a.baseline} vs {a.new})")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
