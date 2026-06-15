#!/usr/bin/env bash
# Run regenerate_tf_dmr.sh for all five lineages, sequentially.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for lineage in mge bach2 epha4 opc ecge; do
  echo "===== $lineage ====="
  bash "$SCRIPT_DIR/regenerate_tf_dmr.sh" "$lineage"
  echo
done

echo "All five lineages complete."