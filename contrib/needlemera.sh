#!/usr/bin/env bash

set -euo pipefail

seq_id="$(pbpaste)"
stumpish revcomp
seq_rc="$(pbpaste)"
stumpish revcomp

needle \
    -asequence asis:"${seq_id}" \
    -bsequence asis:"${seq_rc}" \
    -outfile ~/scratch/needlemera.txt \
    -gapopen 2.0 \
    -gapextend 0.5 \
    -awidth3 1000000


less -S ~/scratch/needlemera.txt

