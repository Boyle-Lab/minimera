minimera_on() {
    "$TESTDIR/../build/minimera" --monotony-fraction 0.001 --output ./out --no-progress --threads 4 "$@"
}

minimera() {
    minimera_on "$TESTDIR/data/example.fastq" "$@"
}

results_on() {
    (head -n1 ./"$1"/foldbacks.csv | cut -d, --complement -f11 && tail -n+2 ./"$1"/foldbacks.csv | cut -d, --complement -f 11 | sort) | column -s, -t
}

results() {
    results_on out "$@"
}
