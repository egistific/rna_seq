#!/usr/bin/env bash
# Creates a local dev environment for this project using micromamba rather
# than plain `conda env create` (the classic conda solver has repeatedly
# gotten OOM-killed on bioconda-channel solves on this machine, regardless
# of available RAM -- micromamba's solver doesn't hit the same issue).
#
# This project doesn't have its own envs/*.yaml files yet (workflow/rules/
# references conda: "../../envs/qc.yaml" etc., but no envs/ directory
# exists)
set -euo pipefail

MAMBA_DIR="${MAMBA_DIR:-.micromamba}"
ENV_NAME="rnaseq"

echo "== 1. Fetch the standalone micromamba binary (no conda install needed) =="
mkdir -p "$MAMBA_DIR"
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xj -C "$MAMBA_DIR" bin/micromamba
export MAMBA_ROOT_PREFIX="$(pwd)/$MAMBA_DIR"
MM="$MAMBA_DIR/bin/micromamba"
"$MM" --version

echo "== 2. Dev env: $ENV_NAME =="
"$MM" create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
  python=3.11 snakemake-minimal=8.4 pandas=2.2 pyyaml "conda>=24.7.1"

cat <<EOF

Done. To use:

  source workflow/scripts/activate_env.sh
EOF
