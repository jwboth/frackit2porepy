#!/usr/bin/env bash
set -euo pipefail

SHARED_DIR="/frackit/shared"
BIN="./frackit2porepy"

echo "[entrypoint] Using shared dir: ${SHARED_DIR}"
echo "[entrypoint] Running: ${BIN}"

# Run binary.
cd frackit/build/appl/frackit2porepy/
"$BIN" "$@"

# Copy known outputs back to shared folder (if they exist)
for f in families.csv disks.csv quads.csv network.brep network.geo; do
  if [[ -f "$f" ]]; then
    echo "[entrypoint] Exporting $f -> ${SHARED_DIR}/$f"
    cp -f "$f" "${SHARED_DIR}/"
  else
    echo "[entrypoint] Note: output file not found: $f"
  fi
done