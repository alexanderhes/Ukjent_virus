# Changelog

What changed for users in each release. Format based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/): MAJOR for changes that break
existing outputs or wrapper calls, MINOR for new features, PATCH for fixes.

Production (`NGS_wrapper.sh`) runs a tagged release (`VERSION="v<version>"`).

## [Unreleased]

## [0.1.1] — 2026-11-19

### Changed
- `--validate_min_contig_len` default 800 → 300. Real viral contigs of
  300–799 bp (gap fillers, low-abundance viruses) were dropped; short artefacts
  are already caught by BLAST and the `verdict_artefact_*` thresholds.
  Contig-based coverage, `blast_assessment` and verdicts can rise; more
  `artefact_suspected` rows appear in the overview.
- `host_moose` profile points to `/mnt/databases/genome_references/moose`.
- `NGS_wrapper.sh` runs `v0.1.1`.

## [0.1.0] — 2026-11-19

First tagged release, not yet feature complete (1.0.0 is kept for that).
The pipeline code equals commit `978654d` on `Assembly_v2` (in production
since 2026-10-12), plus versioning and tooling.

### Added
- Version number (`manifest.version`), printed at the top of every run log and
  in `pipeline_info/report.html`. The wrapper runs the release tag instead of the
  `Assembly_v2` branch and logs the version in the status file.
- Combined genomes per sample and virus (per segment for segmented viruses),
  full length (`<id>_<group>.fasta`, end Ns kept) and trimmed
  (`<id>_<group>_trimmed.fasta`), with contig overhangs beyond the reference
  ends (`--scaffold_max_overhang`, default 300, 0 = off).
- Overview columns `genome_trimmed_len`, `esv_segment_breadth_pct`,
  `blast_ref_cov_bp`, `blast_assessment`; flags `low_rpm`, `low_ref_cov`.
  `genome_assessment` is graded on called reference positions, so overhangs and
  inserted bases no longer raise it.
- Genome plots `<sample>_validation_genomes.pdf` and `<sample>_genome_coverage.pdf`;
  contig plots split into `_validation_contigs.pdf` and `_validation_artefacts.pdf`.
- `--polish_contig_bases` (keep | reads | fallback; default fallback).

### Changed
- **Host filtering is stricter** (`bowtie2 --dovetail -X 1000`): more host
  pairs are removed, so host-removal % rises. EsViritu detections and read
  counts were unchanged in the test samples.
- Consensus: reads decide by majority where depth >= `polish_min_depth` (5);
  below that the contig base is kept (fallback) or N. No IUPAC codes.
- `polish_rounds` default 2 → 1.
- **Breaking:** the genome FASTA under the old file name is now full length
  (end Ns kept); the trimmed genome is `*_trimmed.fasta`.
- Verdict rules: ESV supported needs reads >= 10, RPM >= 10, breadth >= 30 %;
  BLAST confirmed by reference coverage >= 1000 bp or >= 30 % or a partial
  genome. Grades complete / high / partial / low / fragmented.
  Thresholds are still provisional.

### Removed
- `--polish_iupac` and the stats column `n_ambig`.
- `--verdict_contig_min_aln_bp` / `--verdict_contig_min_aln_pct`.

## Before 0.1.0

Untagged. See `git log` up to `978654d`.
