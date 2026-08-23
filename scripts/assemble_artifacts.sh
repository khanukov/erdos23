#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
OUTPUT_DIR="${REPO_ROOT}/.replay-work"
mkdir -p "${OUTPUT_DIR}"

assemble() {
  local target="$1"
  shift
  local parts=("$@")
  if [[ "${#parts[@]}" -eq 0 ]]; then
    echo "no chunks found for ${target}" >&2
    exit 1
  fi
  local temporary
  temporary="$(mktemp "${OUTPUT_DIR}/.assemble.XXXXXX")"
  cat -- "${parts[@]}" > "${temporary}"
  mv -- "${temporary}" "${OUTPUT_DIR}/${target}"
}

assemble "erdos23_global_exact_dual.json" \
  "${REPO_ROOT}"/certificate/erdos23_global_exact_dual.json.part-*
assemble "erdos23_full_replay_bundle_2026-08-23.tar.gz" \
  "${REPO_ROOT}"/release/erdos23_full_replay_bundle_2026-08-23.tar.gz.part-*
assemble "erdos23_full_solution_exact_certificate_2026-08-23.tar.gz" \
  "${REPO_ROOT}"/release/erdos23_full_solution_exact_certificate_2026-08-23.tar.gz.part-*

cd "${REPO_ROOT}"
sha256sum -c ASSEMBLED_SHA256SUMS
echo "ASSEMBLED_ARTIFACTS_OK"
