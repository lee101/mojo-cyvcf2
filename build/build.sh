#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p dist
mojo build --emit shared-lib src/capi.mojo -o dist/libmojo-cyvcf2.so
