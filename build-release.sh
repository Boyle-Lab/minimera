#!/usr/bin/env bash

set -euo pipefail

V="$(git describe --dirty)"


cd build

mkdir -p minimera-"${V}"

tar -czf \
    minimera-"${V}"/minimera-"${V}"-linux-amd64.tar.gz \
    -- \
    minimera minimera.1

cp minimera.sif minimera-"${V}"/minimera-"${V}".sif
