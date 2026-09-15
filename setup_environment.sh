#!/usr/bin/env bash
set -euo pipefail

ENV_FILE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/environment.yml"
ENV_NAME="methylkey"

if command -v mamba >/dev/null 2>&1; then
    CONDA="mamba"
elif command -v conda >/dev/null 2>&1; then
    CONDA="conda"
else
    printf '%s\n' "Error: install Miniforge, Mambaforge, or conda before running this script." >&2
    exit 1
fi

"${CONDA}" env create --file "${ENV_FILE}" --name "${ENV_NAME}"

"${CONDA}" run --name "${ENV_NAME}" Rscript - <<'RSCRIPT'
pak::pkg_install(c(
  "github::pepijn-devries/ggsankeyfier",
  "github::perishky/dmrff",
  "github::sailalithabollepalli/EpiSmokEr",
  "github::achilleasNP/IlluminaHumanMethylationEPICmanifest",
  "github::achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38",
  "github::chiaraherzog/IlluminaMouseMethylationanno.12.v1.mm10",
  "github::IARCbioinfo/methylkey"
))
RSCRIPT

printf '\nEnvironment created: %s\n' "${ENV_NAME}"
printf 'Activate it with: conda activate %s\n' "${ENV_NAME}"