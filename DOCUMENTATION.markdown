# minimera

USAGE: `minimera [OPTIONS] \-\-output=PATH FASTQ`

Minimera is a tool for detecting foldback chimeric and other problematic reads
in Oxford Nanopore data using minimizers.

Reads will be read from the specified fastq file (or stdin if \- is given) and
the results written to foldbacks.csv inside the specified output directory.  If
plots are generated they will be inside a plots/ subdirectory of the output
directory.

The resulting foldbacks.csv file will contain the following columns:

1. read\-id: the ID of the read

2. read\-length: the length (in base pairs) of the read

3. classification: one of foldback/monotonous/normal

4. monotony: monotony score for the read

5. foldback\-point: estimated position of the foldback point for foldbacks, empty
otherwise

6. llqma: length of the longest low\-quality moving average run in the read

7. llqma\-start: start of the longest low\-quality moving average run, blank if no
run exists

8. llqma\-end: end of the longest low\-quality moving average run, blank if no run
exists

9. processing\-time\-microsec: how long minimera took to process this read

New columns may be added after these in future versions of Minimera, so when
processing foldbacks.csv make sure to allow extra columns if you want to remain
forwards compatible.

## Options

*   `-h`, `--help`

    Display help and exit.

*   `-v`, `--version`

    Display version and exit.

*   `-j N`, `--threads=N`

    Number of worker threads to spawn. Does not include the reader and writer
    threads (default: 1).

*   `--progress`

    Report progress as reads are classified (the default).

*   `--no-progress`

    Do not report progress.

*   `-o PATH`, `--output=PATH`

    Output directory (required).


### Foldback Options

*   `-K N`, `--kmer-size=N`

    Size of kmers (in base pairs) to use for minimizers (default: 8).

*   `-W N`, `--window-size=N`

    Size (total) of windows (in base pairs) to use for minimizer sketches
    (default: 16).

*   `-F N`, `--foldback-position-epsilon=N`

    How close (in base pairs) a cluster must be to the beginning of read 2 and
    end of read 1 to be considered as a foldback cluster (default: 50).

*   `-I N`, `--intercept-epsilon=N`

    Epsilon (in y\-intercept space) used to cluster colinear points during the
    initial clustering step (default: 30).

*   `-G N`, `--gap-epsilon=N`

    Maximum width (in y\-intercept space) of an allowable gap in a cluster.  Gaps
    wider than this will result in splitting the cluster (default: 300).

*   `-C N`, `--minimum-cluster-length=N`

    Minimum number of allowable points in a cluster.  Clusters with fewer than
    this many hits will be removed (default: 40).

*   `-A N`, `--minimum-foldback-length-absolute=N`

    Minimum length (in base pairs) of a foldback region.  Regions shorter than
    this will be excluded (default: 50).

*   `-R X`, `--minimum-foldback-length-relative=X`

    Minimum length (as a fraction of total read length) of a foldback region.
    Regions shorter than this will be excluded (default: 0.05).


### Monotony Options

*   `-M X`, `--monotony-threshold=X`

    Monotony score above which reads will be classified as monotonous and not
    analyzed for foldback finding (default: 0.80).


### Low Quality Moving Average Options

*   `-L N`, `--low-quality-window-size=N`

    Size of moving average window when computing LLQMA (default: 10).

*   `-Q Q`, `--low-quality-threshold=Q`

    Threshold (inclusive) for low\-quality Phred scores when computing LLQMA
    (default: 3).


### Plotting Options

Minimera can optionally generate plots of the minimizers for reads.  Generating
these plots is slow, but they can be useful to debug edge cases.

*   `--plot-foldbacks`

    Generate plots for all reads determined to be foldbacks.

*   `--no-plot-foldbacks`

    Do not generate plots for reads determined to be foldbacks (the default).

*   `--plot-normal`

    Generate plots for all reads determined to be normal.

*   `--no-plot-normal`

    Do not generate plots for reads determined to be normal (the default).


## Examples

Run minimera on a FASTQ:

    minimera path/to/data.fastq --output ./results/

Run minimera with 8 worker threads, plotting any foldbacks:

    minimera path/to/data.fastq --threads 8 --plot-foldbacks --output ./results/

