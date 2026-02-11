#!/usr/bin/env bash

set -euo pipefail

# Prep ------------------------------------------------------------------------
version="$(git describe --dirty)"

cd build

# Release directory -----------------------------------------------------------
release_dir="minimera-${version}"

mkdir -p "${release_dir}"
cd "${release_dir}"

# Linux binaries --------------------------------------------------------------
binary_dir="minimera-${version}-linux-amd64"
binary_tar="${binary_dir}.tar.gz"

mkdir -p ./"${binary_dir}"

cp -t ./"${binary_dir}" \
    ../../LICENSE \
    ../minimera \
    ../minimera.1 \
    ../minimera.fish

tar -czf "${binary_tar}" ./"${binary_dir}"

rm -r ./"${binary_dir}"

# Singularity container -------------------------------------------------------
container="minimera-${version}.sif"

cp ../minimera.sif "${container}"
