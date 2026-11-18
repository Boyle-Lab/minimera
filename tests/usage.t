  $ . "$TESTDIR/setup.sh"

Usage:

  $ "$TESTDIR/../build/minimera"
  error: Output directory must be specified.
  [1]

  $ "$TESTDIR/../build/minimera" -h
  minimera - detect foldback chimeric reads (and other artifacts) using minimizers
  
  USAGE: /home/sjl/src/minimera/tests/../build/minimera [OPTIONS] --output=PATH FASTQ
  
  Minimera is a tool for detecting foldback chimeric and other problematic reads
  in Oxford Nanopore data using minimizers.
  
  Reads will be read from the specified fastq file (or stdin if - is given) and
  the results written to foldbacks.csv inside the specified output directory.  If
  plots are generated they will be inside a plots/ subdirectory of the output
  directory.
  
  Options:
    -h, --help            Display help and exit.
    -v, --version         Display version and exit.
    -o PATH, --output PATH
                          Output directory (required).
    -j N, --threads N     Number of worker threads to spawn. Does not include the
                          reader and writer threads (default: 1).
    --progress            Report progress as reads are classified (the default).
    --no-progress         Do not report progress.
  
  Input Options:
    --dna                 Run on DNA input (the default).
    --rna                 Run on RNA input instead of DNA.
    --fastq               Inputs are FASTQ files (the default).
    --bam                 Inputs are BAM files.
  
  Filtering Options:
    Minimera will try to filter out low-quality and/or uninformative reads before
    generating minimizers and performing alignment.
  
    --min-qscore Q        Minimum mean Q-score, reads with a mean Q-score less
                          than this will be classified as failed and not analyzed
                          for foldback finding (default: 9.0).
    --dorado-mean-qscore  Compute mean Q-score using the same method as Dorado
                          (i.e. dropping the first 60 bases first).
    --simple-mean-qscore  Compute mean Q-score as a simple mean of Q-score
                          probabilities (the default).
    --monotony-fraction X Fraction of top minimizers to use when computing
                          monotony score (default: 0.02).
    --monotony-threshold X
                          Monotony score above which reads will be classified as
                          failed and not analyzed for foldback finding (default:
                          0.50).
    --minimum-monotony-length N
                          Minimum read length (in base pairs) for monotony
                          filtering (default: 1000).
  
  Clustering Options:
    For reads that have passed basic quality checks, Minimera will generate
    alignments by clustering minizer matches between the read and its own reverse
    complement.
  
    -k N, --kmer-size N   Size of kmers (in base pairs) to use for minimizers
                          (default: 8).
    -w N, --window-size N Total size of windows (in base pairs) to use for
                          minimizer sketches (default: 16).
    --intercept-epsilon N Epsilon (in y-intercept space) used to cluster colinear
                          points during the initial clustering step (default: 30).
    --gap-epsilon N       Maximum width (in y-intercept space) of an allowable gap
                          in a cluster.  Gaps wider than this will result in
                          splitting the cluster (default: 300).
    --cluster-lis-fraction X
                          Fraction of points that must be included in the LIS of a
                          cluster for it to be valid (default: 0.10).
    --minimum-cluster-length N
                          Minimum number of allowable points in a cluster.
                          Clusters with fewer than this many hits will be removed
                          (default: 40).
  
  Foldback Options:
    Minimera will use the alignment clusters to classify reads as normal or
    foldback chimeric reads.
  
    --foldback-position-epsilon N
                          How close (in base pairs) a cluster must be to the
                          beginning/end of the reads to be considered as a
                          foldback cluster (default: 80).
    --minimum-foldback-length-absolute N
                          Minimum length (in base pairs) of a foldback region.
                          Regions shorter than this will be excluded (default:
                          50).
    --minimum-foldback-length-relative X
                          Minimum length (as a fraction of total read length) of a
                          foldback region.  Regions shorter than this will be
                          excluded (default: 0.05).
  
  Low Quality Region Options:
    Minimera computes the longest low-quality region (if any) of each read.
  
    --low-quality-window-size N
                          Size of moving average window when computing LLQR
                          (default: 10).
    --low-quality-threshold Q
                          Threshold (inclusive) for low-quality Q-scores when
                          computing LLQR (default: 3).
  
  Plotting Options:
    Minimera can optionally generate plots of the minimizers for reads.
    Generating these plots is slow, but they can be useful to debug edge cases.
  
    --plot-foldbacks      Generate plots for all reads classified as foldbacks.
    --no-plot-foldbacks   Do not generate plots for reads classified as foldbacks
                          (the default).
    --plot-normal         Generate plots for all non-foldback reads.
    --no-plot-normal      Do not generate plots for non-foldback reads (the
                          default).
  
  Examples:
  
    Run minimera on a FASTQ:
  
        minimera path/to/data.fastq --output ./results/
  
    Run minimera with 8 worker threads, plotting any foldbacks:
  
        minimera path/to/data.fastq --threads 8 --plot-foldbacks --output ./results/
