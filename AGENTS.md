# Repository Guidelines

## Project Structure & Module Organization

The command-line entry point is `script/dbitm.sh`. Pipeline stages live in
`script/steps/`: numbered shell wrappers coordinate tools, while supporting
Python implementations are under `script/steps/python/`. Keep stage numbers and
names aligned, for example `07.saturation.sh` and `python/07.saturation.py`.

Configuration defaults belong in `config/dbitm.config.example.sh`; the local
`config/dbitm.config.sh` is ignored and may contain machine-specific paths.
Barcode manifests and technical documentation are under `docs/`. Exploratory
or project-specific downstream work belongs in `analysis/`, not in pipeline
stage code. This repository currently has no dedicated automated test suite.

## Build, Test, and Development Commands

Use Pixi for all environments and commands:

```bash
pixi install
pixi run init
pixi run -e default bash -n script/steps/07.saturation.sh
pixi run -e default python script/steps/python/07.saturation.py --help
dbitm smc all --input /path/to/sample/fastq --dry-run
```

`pixi install` resolves the locked environment. `pixi run init` installs the
`dbitm` launcher. Shell syntax checks and `--help` catch lightweight regressions;
the pipeline dry-run validates configuration and stage dependencies without
writing results. Check `pixi.toml` before assuming additional tasks exist.

## Coding Style & Naming Conventions

Use four-space indentation, type hints, `snake_case` functions and variables,
and `UPPER_CASE` constants in Python. Prefer `pathlib.Path`, explicit validation,
and atomic output replacement for generated files. Shell scripts should use
Bash, `set -euo pipefail`, quoted variables, arrays for command arguments, and
clear `[dbitm]` log prefixes. Avoid unrelated formatting or refactoring.

## Testing Guidelines

Design a small test around the changed behavior. Run Python and bioinformatics
tools through `pixi run -e default`. Put synthetic FASTQ, BAM, and coverage
fixtures under an isolated `/tmp` directory and remove them after inspection.
For workflow changes, test the affected assay and stage with `--dry-run`; for
plot changes, inspect the generated image and its summary TSV together.

## Commit & Pull Request Guidelines

Follow the repository's Conventional Commit history, for example
`fix(qc): align SmC read inputs` or `feat(smc): add strand-aware filtering`.
Do not commit generated results, local configuration, or reference paths.

Pull requests should describe the affected assays and stages, explain behavior
and configuration changes, list checks performed, and link relevant issues.
Include before/after plots when QC visualization changes. Keep commits focused
and call out anything not tested on full data.

## Security & Data Handling

Never commit credentials, patient data, raw sequencing files, reference genomes,
or absolute infrastructure paths. Preserve existing outputs and uncommitted
work; use scratch or `/tmp` for disposable intermediates.
