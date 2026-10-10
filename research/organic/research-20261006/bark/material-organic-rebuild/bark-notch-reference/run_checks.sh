#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
ulimit -v 1048576
timeout 115s python tests/test_reference.py > receipts/tests.log 2>&1
timeout 115s python source/run_reference.py > receipts/reference.log 2>&1
timeout 115s python source/verify_receipts.py
cat receipts/tests.log receipts/reference.log
