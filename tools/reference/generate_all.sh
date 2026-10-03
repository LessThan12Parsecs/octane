#!/usr/bin/env bash
# Regenerate every oracle fixture.
set -euo pipefail
cd "$(dirname "$0")/../.."
for f in tools/reference/*.py; do
  echo "== $f"
  .venv/bin/python "$f"
done
