#!/usr/bin/env bash

set -euo pipefail

function fastq {
    echo "@scratch"
    pbpaste | head -n1
    echo
    echo "+"
    pbpaste | head -n1 | sed -Ee 's/./I/g'
    echo
}


scratch="/home/sjl/scratch/minimera-clipboard-scratch"
mkdir -p "$scratch"

if test -f "$scratch"/plots/scratch.png; then
    rm "$scratch"/plots/scratch.png
fi

minimera <(fastq) \
    --output "$scratch" \
    --plot-foldbacks \
    --plot-normal

qimgv "$scratch"/plots/scratch.png
