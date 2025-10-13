#!/usr/bin/env bash

set -eu

set -x
function seq_to_fastq {
    echo "@scratch"
    pbpaste | tr -d '\n' | head -n1
    echo
    echo "+"
    pbpaste | tr -d '\n' | head -n1 | sed -Ee 's/./I/g'
    echo
}

function clean_fastq {
    pbpaste | awk '
    NR <= 4 {
        printf("%s\n", $0)
    }
    '
}

function read_id {
    pbpaste | head -n1 | cut -f1 | tail -c +2
}

scratch="/home/sjl/scratch/minimera-clipboard-scratch"
mkdir -p "$scratch"

# if test -f "$scratch"/plots/scratch.png; then
#     rm "$scratch"/plots/scratch.png
# fi

id=$(read_id)
# fastq

minimera <(clean_fastq) \
    --output "$scratch" \
    --plot-foldbacks \
    --plot-normal

cp  "$scratch"/plots/"$id".png \
    "$scratch"/plots/scratch.png

if ps -ax | grep -P '[q]imgv' > /dev/null; then
    echo qimgv already running
else
    nohup qimgv "$scratch"/plots/scratch.png >/dev/null &
fi
