#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

sha256sum -c ARTIFACT_SHA256SUMS

./scripts/assemble_artifacts.sh

COMPACT_ARCHIVE=".replay-work/erdos23_full_solution_exact_certificate_2026-08-23.tar.gz"
FULL_ARCHIVE=".replay-work/erdos23_full_replay_bundle_2026-08-23.tar.gz"

gzip -t "${COMPACT_ARCHIVE}"

ARCHIVE_LIST="$(tar -tzf "${FULL_ARCHIVE}")"
ROOT_COUNT="$(printf '%s\n' "${ARCHIVE_LIST}" | awk '/\/root_[0-9][0-9][0-9]\.npz$/ {count++} END {print count+0}')"
if [[ "${ROOT_COUNT}" != "107" ]]; then
  echo "expected 107 K7 root caches, found ${ROOT_COUNT}" >&2
  exit 1
fi

TEMP_COUNT="$(printf '%s\n' "${ARCHIVE_LIST}" | awk '/\/\.root_/ {count++} END {print count+0}')"
if [[ "${TEMP_COUNT}" != "0" ]]; then
  echo "unexpected temporary K7 cache files: ${TEMP_COUNT}" >&2
  exit 1
fi

python3 -m py_compile verifier/*.py generation/*.py

echo "ARTIFACT_INTEGRITY_OK"
echo "k7_root_caches=${ROOT_COUNT}"
