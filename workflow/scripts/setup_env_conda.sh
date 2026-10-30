#!/usr/bin/env bash
# Creates a local dev environment for this project using plain `conda`
# instead of micromamba (see setup_env.sh for the micromamba version).
#
# The classic conda solver has previously gotten OOM-killed on
# bioconda-channel solves on this machine -- --solver=libmamba avoids that.
# conda-libmamba-solver is already installed in this machine's base conda,
# so no extra setup is needed for the flag to work.
#
# This project doesn't have its own envs/*.yaml files yet (workflow/rules/
# references conda: "../../envs/qc.yaml" etc., but no envs/ directory
# exists)
set -euo pipefail

ENV_NAME="rnaseq"

echo "== Dev env: $ENV_NAME =="
conda create -y -n "$ENV_NAME" --solver=libmamba -c conda-forge -c bioconda \
  python=3.11 snakemake-minimal=8.4 pandas=2.2 pyyaml "conda>=24.7.1"

cat <<EOF

Done. To use:

  source workflow/scripts/activate_env_conda.sh
EOF
