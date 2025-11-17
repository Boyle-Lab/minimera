#!/usr/bin/env bash

set -euo pipefail

V=$(git describe --dirty)

buildapp \
        --load-system 'minimera' \
        --eval '(setf minimera::*version* "'"${V}"'")' \
        --entry 'minimera:toplevel' \
        --manifest-file 'build/asdf-manifest' \
        --compress-core \
        --output 'build/minimera'
