#!/usr/bin/env python3
"""
Collects read-count statistics from the fastp JSON outputs and HOST_FILTER.

Pipeline order: FASTP_TRIM -> FASTP_DEDUP -> HOST_FILTER, so
    raw_reads           = FASTP_TRIM input
    trimmed_reads       = FASTP_TRIM output  (= FASTP_DEDUP input)
    dedup_reads         = FASTP_DEDUP output (= HOST_FILTER input)
    host_filtered_reads = HOST_FILTER output (= reads analysed by EsViritu/SPAdes)
All counts are R1 + R2 reads.

Usage:
    collect_read_stats.py <sample_id> <fastp_trim_json> <fastp_dedup_json> <hostfilt_count_txt> <out_tsv>

Output columns:
    sample_ID | raw_reads | trimmed_reads | dedup_reads | host_filtered_reads
    trim_removed_pct | dup_rate_pct | host_removal_pct
    (each percentage is relative to the step's input: raw, trimmed and dedup reads)
"""

import json
import sys

sample_id, fastp_trim_json, fastp_dedup_json, hostfilt_count_txt, out_tsv = sys.argv[1:]

with open(fastp_trim_json) as f:
    trim = json.load(f)

with open(fastp_dedup_json) as f:
    dedup = json.load(f)

with open(hostfilt_count_txt) as f:
    host_filtered = int(f.read().strip() or 0)

raw_reads   = trim['summary']['before_filtering']['total_reads']
trimmed     = trim['summary']['after_filtering']['total_reads']
dedup_reads = dedup['summary']['after_filtering']['total_reads']

trim_removed_pct = round((raw_reads - trimmed)         / raw_reads * 100, 2) if raw_reads > 0 else 0.0
dup_rate_pct     = round((trimmed - dedup_reads)       / trimmed * 100, 2) if trimmed > 0 else 0.0
host_removal_pct = round((dedup_reads - host_filtered) / dedup_reads * 100, 2) if dedup_reads > 0 else 0.0

header = "\t".join([
    "sample_ID", "raw_reads", "trimmed_reads", "dedup_reads", "host_filtered_reads",
    "trim_removed_pct", "dup_rate_pct", "host_removal_pct"
])
row = "\t".join(str(x) for x in [
    sample_id, raw_reads, trimmed, dedup_reads, host_filtered,
    trim_removed_pct, dup_rate_pct, host_removal_pct
])

with open(out_tsv, 'w') as f:
    f.write(header + "\n" + row + "\n")
