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
| `-d` | EsViritu database alias | `v3_2_4` |
| `--resume` | Resume a previous Nextflow run | off |

### Available database aliases
- `v3_2_4` — full EsViritu database (viral metagenomics, SPAdes `--meta`)
- `HEV` — HEV-specific database (SPAdes `--rnaviral` selected automatically)

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
```

# Follow the live wrapper log from outside the screen
tail -f /home/ngs/esv_wrapper.log

# Check the last status of a specific run
cat ~/esv_test_run_status.txt
```
