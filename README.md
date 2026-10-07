# EsViritu Pipeline — Running the Wrapper

### Required arguments

| Flag | Description | Example |
|------|-------------|---------|
| `-r` | Run name | `NGS_SEQ-20260210-01` |
| `-a` | Agens subfolder on the N-drive | `UkjentVirus` |
| `-y` | Year | `2026` |

### Optional arguments

| Flag | Description | Default |
|------|-------------|---------|
| `-d`, `--db` | EsViritu database alias | `v3_2_4` |
| `-H`, `--host` | Host reference alias used for host read removal | `human` |
| `--sensitive-filter` | High-sensitivity host filtering (slower; R1 and R2 mapped separately with bowtie2 `--very-sensitive-local`, then re-paired). Default is a single paired-end bowtie2 `--sensitive-local` run. | off |
| `--resume` | Resume a previous Nextflow run | off |
| `-h`, `--help` | Show usage | |

### Available database aliases
- `v3_2_4` — full EsViritu database (viral metagenomics, SPAdes `--meta`)
- `HEV` — HEV-specific database (SPAdes `--rnaviral` selected automatically)

### Available host aliases
- `human` — human T2T (GCA_009914755.4) + PhiX174
- `moose` — Norwegian moose (*Alces alces*) — index not built yet; the profile path in `nextflow.config` is a placeholder

New database or host aliases need both a `db_<alias>` / `host_<alias>` profile in
`nextflow.config` and an entry in `VALID_DB_ALIASES` / `VALID_HOST_ALIASES` in `NGS_wrapper.sh`.

## Examples

```bash
# Standard metagenomics run
screen -S Test_run -d -m bash /home/ngs/ngs_scripts/ukjent_virus/NGS_wrapper.sh \
-r test_run \
-a UkjentVirus \
-y 2026

# HEV-specific run
screen -S Test_run -d -m bash /home/ngs/ngs_scripts/ukjent_virus/NGS_wrapper.sh \
-r test_run \
-a UkjentVirus \
-y 2026 \
-d HEV

# High-sensitivity host filtering (slower)
screen -S Test_run -d -m bash /home/ngs/ngs_scripts/ukjent_virus/NGS_wrapper.sh \
-r test_run \
-a UkjentVirus \
-y 2026 \
--sensitive-filter
```

## Monitoring

```bash
# Follow the live wrapper log from outside the screen
tail -f /home/ngs/esv_wrapper.log

# Check the last status of a specific run (esv_<run name>_status.txt)
cat ~/esv_test_run_status.txt
```
