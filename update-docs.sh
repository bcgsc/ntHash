#!/usr/bin/env bash

if [ -z "${MESON_SOURCE_ROOT}" ]; then
  echo "[ERROR] This script can only be ran with meson!"
  exit 1
fi

set -euo pipefail

cd "${MESON_SOURCE_ROOT}"
rm -r docs/
doxygen