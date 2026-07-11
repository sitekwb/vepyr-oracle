#!/usr/bin/env bash
# Sync the oracle package to the cluster harness dir.
set -euo pipefail
cd "$(dirname "$0")"
rsync -av --delete \
  --exclude '__pycache__' --exclude '.pytest_cache' --exclude 'legacy' --exclude 'tests' \
  oracle bin slurm \
  hpc:vepyr/work/
echo "[deploy] synced oracle/ bin/ slurm/ -> hpc:~/vepyr/work/"
