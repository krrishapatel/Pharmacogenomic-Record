#!/usr/bin/env bash
set -euo pipefail
VERSION="3.4.0"
DEST="data/pharmcat_positions_${VERSION}.vcf"
URL="https://github.com/PharmGKB/PharmCAT/releases/download/v${VERSION}/pharmcat_positions_${VERSION}.vcf"
mkdir -p data
curl -sSL -o "$DEST" "$URL"
echo "wrote $DEST ($(wc -c < "$DEST") bytes)"
