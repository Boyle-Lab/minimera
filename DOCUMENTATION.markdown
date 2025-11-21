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

6. llq: length of the longest low\-quality region in the read

7. llq\-start: start of the longest low\-quality region, blank if none exists

8. llq\-end: end of the longest low\-quality region, blank if none exists

9. processing\-time\-microsec: how long minimera took to process this read

New columns may be added after these in future versions of Minimera, so when
processing foldbacks.csv make sure to allow extra columns if you want to remain
forwards compatible.

## Options

*   `-h`, `--help`

    Display help and exit.

*   `-v`, `--version`

    Display version and exit.

*   `-o PATH`, `--output=PATH`

    Output directory (required).

*   `-j N`, `--threads=N`

    Number of worker threads to spawn. Does not include the reader and writer
    threads (default: 1).

*   `--progress`

    Report progress as reads are classified (the default).

*   `--no-progress`

    Do not report progress.


### Filtering Options

Minimera will try to filter out low\-quality and/or uninformative reads before
generating minimizers.

*   `--min-qscore=Q`

    Minimum mean Q\-score, reads with a mean Q\-score less than this will be
    classified as failed and not analyzed for foldback finding (default: 9.0).

*   `--monotony-threshold=X`

    Monotony score above which reads will be classified as failed and not
    analyzed for foldback finding (default: 0.50).


### Foldback Options

For reads that have passed basic quality checks, Minimera will classify them as
normal or foldback chimeric reads.

*   `-K N`, `--kmer-size=N`

    Size of kmers (in base pairs) to use for minimizers (default: 8).

*   `-W N`, `--window-size=N`

    Total size of windows (in base pairs) to use for minimizer sketches
    (default: 16).

*   `--foldback-position-epsilon=N`

    How close (in base pairs) a cluster must be to the beginning of read 2 and
    end of read 1 to be considered as a foldback cluster (default: 50).

*   `--intercept-epsilon=N`

    Epsilon (in y\-intercept space) used to cluster colinear points during the
    initial clustering step (default: 30).

*   `--gap-epsilon=N`

    Maximum width (in y\-intercept space) of an allowable gap in a cluster.  Gaps
    wider than this will result in splitting the cluster (default: 300).

*   `--minimum-cluster-length=N`

    Minimum number of allowable points in a cluster.  Clusters with fewer than
    this many hits will be removed (default: 40).

*   `--minimum-foldback-length-absolute=N`

    Minimum length (in base pairs) of a foldback region.  Regions shorter than
    this will be excluded (default: 50).

*   `--minimum-foldback-length-relative=X`

    Minimum length (as a fraction of total read length) of a foldback region.
    Regions shorter than this will be excluded (default: 0.05).


### Low Quality Region Options

Minimera computes the longest low\-quality region (if any) of each read.

*   `--low-quality-window-size=N`

    Size of moving average window when computing LLQR (default: 10).

*   `--low-quality-threshold=Q`

    Threshold (inclusive) for low\-quality Q\-scores when computing LLQR (default:
    3).


### Plotting Options

Minimera can optionally generate plots of the minimizers for reads.  Generating
these plots is slow, but they can be useful to debug edge cases.

*   `--plot-foldbacks`

    Generate plots for all reads classified as foldbacks.

*   `--no-plot-foldbacks`

    Do not generate plots for reads classified as foldbacks (the default).

*   `--plot-normal`

    Generate plots for all non\-foldback reads.

*   `--no-plot-normal`

    Do not generate plots for non\-foldback reads (the default).


## Examples

Run minimera on a FASTQ:

    minimera path/to/data.fastq --output ./results/

Run minimera with 8 worker threads, plotting any foldbacks:

    minimera path/to/data.fastq --threads 8 --plot-foldbacks --output ./results/

