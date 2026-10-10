#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
ulimit -v 1048576
timeout 115s python source/run_tied_refinement.py --protocol ADDITIONAL_REFINEMENT_64.json > receipts/mesh_64_tied.log 2>&1
cat receipts/mesh_64_tied.log
