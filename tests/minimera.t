  $ minimera() { "$TESTDIR/../build/minimera" "$TESTDIR/data/example.fastq" --output ./out --no-progress "$@"; }
  $ results() { (head -n1 ./out/foldbacks.csv | cut -d, --complement -f11 && tail -n+2 ./out/foldbacks.csv | cut -d, --complement -f 11 | sort) | column -s, -t ; }

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
  
  Filtering Options:
    Minimera will try to filter out low-quality and/or uninformative reads before
    generating minimizers.
  
    --min-qscore Q        Minimum mean Q-score, reads with a mean Q-score less
                          than this will be classified as failed and not analyzed
                          for foldback finding (default: 9.0).
    --dorado-mean-qscore  Compute mean Q-score using the same method as Dorado
                          (i.e. dropping the first 60 bases first).
    --simple-mean-qscore  Compute mean Q-score as a simple mean of Q-score
                          probabilities (the default).
    --monotony-threshold X
                          Monotony score above which reads will be classified as
                          failed and not analyzed for foldback finding (default:
                          0.50).
  
  Foldback Options:
    For reads that have passed basic quality checks, Minimera will classify them
    as normal or foldback chimeric reads.
  
    -k N, --kmer-size N   Size of kmers (in base pairs) to use for minimizers
                          (default: 8).
    -w N, --window-size N Total size of windows (in base pairs) to use for
                          minimizer sketches (default: 16).
    --foldback-position-epsilon N
                          How close (in base pairs) a cluster must be to the
                          beginning of read 2 and end of read 1 to be considered
                          as a foldback cluster (default: 80).
    --intercept-epsilon N Epsilon (in y-intercept space) used to cluster colinear
                          points during the initial clustering step (default: 30).
    --gap-epsilon N       Maximum width (in y-intercept space) of an allowable gap
                          in a cluster.  Gaps wider than this will result in
                          splitting the cluster (default: 300).
    --minimum-cluster-length N
                          Minimum number of allowable points in a cluster.
                          Clusters with fewer than this many hits will be removed
                          (default: 40).
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

Run on a FASTQ with some examples:

  $ minimera
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore         llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.324421040730403  0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.346002226031022  0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.358222597036377  0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.344128583944823  0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.368226451679934  0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.393533364827064  0                           0
  foldback-longboi                      127949       foldback        0.0069    63970           19.3547015335795    0                           0
  foldback-lqr-end                      4320         foldback        0.0051    1996            11.580066231058375  323   3997        4320      323
  foldback-lqr-mid-gap                  4620         foldback        0.0051    1993            9.33420463893992    626   1997        2623      626
  foldback-lqr-post-gap                 4229         foldback        0.0051    1993            12.683751180511019  235   3497        3732      235
  foldback-lqr-pre-gap                  4229         foldback        0.0051    2228            12.683751180510981  235   497         732       235
  foldback-lqr-start                    4320         foldback        0.0051    2316            11.580066231058726  323   0           323       323
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.387276545523914  0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.31231293285616   0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.32442732791359   0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.365481608591153  0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.391866552150862  0                           0
  longboi                               120000       normal          0.0067                    19.35712473520419   0                           0
  low-quality-regions                   712          normal          0.0168                    9.501463163531792   102   443         545       122
  monotonous-read                       60007        failed          0.7926                    19.524198987162304  0                           0
  non-foldback                          2800         normal          0.0073                    19.355641445904553  0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.350417940453735  0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.344804016260536  0                           0

Disable LQR computation and splicing by making the threshold 0:

  $ minimera --low-quality-threshold 0
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore         llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.324421040730403  0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.346002226031022  0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.358222597036377  0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.344128583944823  0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.368226451679934  0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.393533364827064  0                           0
  foldback-longboi                      127949       foldback        0.0069    63970           19.3547015335795    0                           0
  foldback-lqr-end                      4320         normal          0.0047                    11.580066231058375  0                           0
  foldback-lqr-mid-gap                  4620         normal          0.0044                    9.33420463893992    0                           0
  foldback-lqr-post-gap                 4229         normal          0.0048                    12.683751180511019  0                           0
  foldback-lqr-pre-gap                  4229         normal          0.0048                    12.683751180510981  0                           0
  foldback-lqr-start                    4320         foldback        0.0047    2316            11.580066231058726  0                           0
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.387276545523914  0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.31231293285616   0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.32442732791359   0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.365481608591153  0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.391866552150862  0                           0
  longboi                               120000       normal          0.0067                    19.35712473520419   0                           0
  low-quality-regions                   712          normal          0.0208                    9.501463163531792   0                           0
  monotonous-read                       60007        failed          0.7926                    19.524198987162304  0                           0
  non-foldback                          2800         normal          0.0073                    19.355641445904553  0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.350417940453735  0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.344804016260536  0                           0
