#!/usr/bin/env python3
"""
Collects read-count statistics from HOST_FILTER and fastp JSON outputs.

Usage:
    collect_read_stats.py <sample_id> <raw_reads_txt> <fastp_trim_json> <fastp_dedup_json> <out_tsv>

Output columns:
    sample_ID | raw_reads | host_filtered_reads | trimmed_reads | dedup_reads
    host_removal_pct | trim_removed_pct | dup_rate_pct
"""

import json
import sys

sample_id, raw_reads_txt, fastp_trim_json, fastp_dedup_json, out_tsv = sys.argv[1:]

with open(raw_reads_txt) as f:
    raw_reads = int(f.read().strip())

with open(fastp_trim_json) as f:
    trim = json.load(f)

with open(fastp_dedup_json) as f:
    dedup = json.load(f)

host_filtered = trim['summary']['before_filtering']['total_reads']
trimmed       = trim['summary']['after_filtering']['total_reads']
dedup_reads   = dedup['summary']['after_filtering']['total_reads']

host_removal_pct = round((raw_reads - host_filtered) / raw_reads * 100, 2) if raw_reads > 0 else 0.0
trim_removed_pct = round((host_filtered - trimmed)   / host_filtered * 100, 2) if host_filtered > 0 else 0.0
dup_rate_pct     = round((trimmed - dedup_reads)      / trimmed * 100, 2) if trimmed > 0 else 0.0

header = "\t".join([
    "sample_ID", "raw_reads", "host_filtered_reads",
    "trimmed_reads", "dedup_reads",
    "host_removal_pct", "trim_removed_pct", "dup_rate_pct"
])
row = "\t".join(str(x) for x in [
    sample_id, raw_reads, host_filtered,
    trimmed, dedup_reads,
    host_removal_pct, trim_removed_pct, dup_rate_pct
])

with open(out_tsv, 'w') as f:
    f.write(header + "\n" + row + "\n")
