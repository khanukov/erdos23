#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
ARCHIVE="${REPO_ROOT}/.replay-work/erdos23_full_replay_bundle_2026-08-23.tar.gz"
PYTHON_BIN="${PYTHON_BIN:-python3}"
REPLAY_JOBS="${REPLAY_JOBS:-4}"

"${REPO_ROOT}/scripts/assemble_artifacts.sh"

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/erdos23-replay.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT

tar -xzf "${ARCHIVE}" -C "${WORK_DIR}"
BUNDLE_ROOT="${WORK_DIR}/erdos23_full_replay_bundle"

cd "${BUNDLE_ROOT}"
sha256sum -c BUNDLE_SHA256SUMS

"${PYTHON_BIN}" verify_moment_vector.py \
  --flagsdp ./public_flagsdp_data \
  --public-anc ./public_anc \
  --jobs "${REPLAY_JOBS}"

for CERTIFICATE in erdos23_global_exact_dual.json erdos23_global_exact_dual.pkl; do
  "${PYTHON_BIN}" verify_exact_fixed_gram_dual.py \
    --certificate "./${CERTIFICATE}" \
    --flagsdp ./public_flagsdp_data \
    --public-anc ./public_anc \
    --k7-cache ./k7_compact_v1 \
    --target-n 1
done

echo "FULL_REPLAY_OK"
