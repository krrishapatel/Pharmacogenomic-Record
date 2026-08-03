#!/usr/bin/env bash
set -euo pipefail
# The reference table is package data, so it lives inside the package tree and
# ships in the wheel -- not in a top-level data/ dir the installed tool could
# never see. Keep this path in step with tool.setuptools.package-data.
VERSION="3.4.0"
EXPECTED_BYTES=64934
DATA_DIR="src/pharmacogenomic_record/data"
DEST="${DATA_DIR}/pharmcat_positions_${VERSION}.vcf"
URL="https://github.com/PharmGKB/PharmCAT/releases/download/v${VERSION}/pharmcat_positions_${VERSION}.vcf"
mkdir -p "$DATA_DIR"

# Download to a temp file and validate before replacing the pinned reference.
# Without -f, curl writes the HTTP error body to the output file and exits 0,
# which would silently overwrite the good file with "Not Found".
TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT
curl -fsSL -o "$TMP" "$URL"

ACTUAL_BYTES="$(wc -c < "$TMP" | tr -d " ")"
if [ "$ACTUAL_BYTES" -ne "$EXPECTED_BYTES" ]; then
    echo "refusing to install: expected ${EXPECTED_BYTES} bytes, got ${ACTUAL_BYTES}" >&2
    exit 1
fi

mv "$TMP" "$DEST"
trap - EXIT
echo "wrote $DEST (${ACTUAL_BYTES} bytes)"
