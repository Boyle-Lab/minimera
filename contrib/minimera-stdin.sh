#!/usr/bin/env bash

# Copyright 2025 Steve Losh and contributors
# SPDX-License-Identifier: GPL-3.0-or-later

set -eu

set -x

input=$(cat)

function clean_fastq {
    echo "${input}" | awk '
    NR <= 4 {
        printf("%s\n", $0)
    }
    '
}

function read_id {
    echo "${input}" | head -n1 | cut -d ' ' -f1 | tail -c +2
}

scratch="/home/sjl/scratch/minimera-scratch"
mkdir -p "$scratch"

id=$(read_id)

minimera <(clean_fastq) \
    --output "$scratch" \
    --plot-foldbacks \
    --plot-normal \
    --min-qscore 0 \
    "$@"


cp  "$scratch"/plots/"$id".png \
    "$scratch"/plots/scratch.png

if ps -ax | grep -P '[q]imgv' > /dev/null; then
    echo qimgv already running
else
    nohup qimgv "$scratch"/plots/scratch.png >/dev/null &
fi
