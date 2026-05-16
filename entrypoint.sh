#!/usr/bin/env bash
set -euo pipefail

SHARED_DIR="/frackit/shared"
BIN="./frackit2porepy"

echo "[entrypoint] Using shared dir: ${SHARED_DIR}"
echo "[entrypoint] Running: ${BIN}"

# Pass an input file from $SHARED_DIR.
INPUT="${SHARED_DIR}/config.toml"
# exec "$BIN" --input "$INPUT"

# Run binary.
cd frackit/build/appl/frackit2porepy/
"$BIN" "$@"

# Copy known outputs back to shared folder (if they exist)
for f in disks.csv network.brep network.geo; do
  if [[ -f "$f" ]]; then
    echo "[entrypoint] Exporting $f -> ${SHARED_DIR}/$f"
    cp -f "$f" "${SHARED_DIR}/"
  else
    echo "[entrypoint] Note: output file not found: $f"
  fi
done