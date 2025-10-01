#!/usr/bin/env bash

set -euo pipefail

function fastq {
    echo "@00000000-0000-0000-0000-000000000000"
    pbpaste | head -n1
    echo
    echo "+"
    pbpaste | head -n1 | sed -Ee 's/./I/g'
    echo
}


scratch="/home/sjl/scratch/minimera-clipboard-scratch"
mkdir -p "$scratch"

minimera <(fastq) \
    --output "$scratch" \
    --plot-foldbacks \
    --plot-normal

qimgv "$scratch"/plots/00000000-0000-0000-0000-000000000000.png
