  $ minimera_on() { "$TESTDIR/../build/minimera" --monotony-fraction 0.001 --output ./out --no-progress --threads 4 "$@"; }
  $ minimera() { minimera_on "$TESTDIR/data/example.fastq" "$@"; }
  $ results_on() { (head -n1 ./"$1"/foldbacks.csv | cut -d, --complement -f11 && tail -n+2 ./"$1"/foldbacks.csv | cut -d, --complement -f 11 | sort) | column -s, -t ; }
  $ results() { results_on out "$@"; }

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
    generating minimizers.
  
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
  
  Foldback Options:
    For reads that have passed basic quality checks, Minimera will classify them
    as normal or foldback chimeric reads.
  
    -k N, --kmer-size N   Size of kmers (in base pairs) to use for minimizers
                          (default: 8).
    -w N, --window-size N Total size of windows (in base pairs) to use for
                          minimizer sketches (default: 16).
    --foldback-position-epsilon N
                          How close (in base pairs) a cluster must be to the
                          beginning/end of the reads to be considered as a
                          foldback cluster (default: 80).
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
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.3244      0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.3460      0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.3582      0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.3441      0                           0
  big-deletion                          3900         normal          0.0052                    19.3885      0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.3682      0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.3935      0                           0
  foldback-longboi                      127949       foldback        0.0072    63970           19.3547      0                           0
  foldback-lqr-end                      4320         foldback        0.0051    1996            11.5801      323   3997        4320      323
  foldback-lqr-mid-gap                  4620         foldback        0.0051    1993            9.3342       626   1997        2623      626
  foldback-lqr-post-gap                 4229         foldback        0.0051    1993            12.6838      235   3497        3732      235
  foldback-lqr-pre-gap                  4229         foldback        0.0051    2228            12.6838      235   497         732       235
  foldback-lqr-start                    4320         foldback        0.0051    2316            11.5801      323   0           323       323
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.3873      0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.3123      0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.3244      0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.3655      0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.3919      0                           0
  large-end-gap-foldback                1720         foldback        0.0119    796             19.3850      0                           0
  longboi                               120000       normal          0.0067                    19.3571      0                           0
  low-quality-first-bases               860          normal          0.0242                    20.5177      0                           0
  low-quality-regions                   712          normal          0.0168                    9.5015       102   443         545       122
  monotonous-read                       60007        failed          0.7926                    19.5242      0                           0
  non-foldback                          2800         normal          0.0073                    19.3556      0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.3504      0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.3448      0                           0
  small-foldback-160bp-total            1160         normal          0.0178                    19.3520      0                           0

Disable LQR computation and splicing by making the threshold 0:

  $ minimera --low-quality-threshold 0
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.3244      0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.3460      0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.3582      0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.3441      0                           0
  big-deletion                          3900         normal          0.0052                    19.3885      0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.3682      0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.3935      0                           0
  foldback-longboi                      127949       foldback        0.0072    63970           19.3547      0                           0
  foldback-lqr-end                      4320         foldback        0.0047    1996            11.5801      0                           0
  foldback-lqr-mid-gap                  4620         normal          0.0044                    9.3342       0                           0
  foldback-lqr-post-gap                 4229         normal          0.0048                    12.6838      0                           0
  foldback-lqr-pre-gap                  4229         normal          0.0048                    12.6838      0                           0
  foldback-lqr-start                    4320         foldback        0.0047    2316            11.5801      0                           0
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.3873      0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.3123      0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.3244      0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.3655      0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.3919      0                           0
  large-end-gap-foldback                1720         foldback        0.0119    796             19.3850      0                           0
  longboi                               120000       normal          0.0067                    19.3571      0                           0
  low-quality-first-bases               860          normal          0.0242                    20.5177      0                           0
  low-quality-regions                   712          normal          0.0208                    9.5015       0                           0
  monotonous-read                       60007        failed          0.7926                    19.5242      0                           0
  non-foldback                          2800         normal          0.0073                    19.3556      0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.3504      0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.3448      0                           0
  small-foldback-160bp-total            1160         normal          0.0178                    19.3520      0                           0

Mean qscore filtering:

  $ minimera --min-qscore 10
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.3244      0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.3460      0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.3582      0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.3441      0                           0
  big-deletion                          3900         normal          0.0052                    19.3885      0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.3682      0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.3935      0                           0
  foldback-longboi                      127949       foldback        0.0072    63970           19.3547      0                           0
  foldback-lqr-end                      4320         foldback        0.0051    1996            11.5801      323   3997        4320      323
  foldback-lqr-mid-gap                  4620         failed          0.0051                    9.3342       626   1997        2623      626
  foldback-lqr-post-gap                 4229         foldback        0.0051    1993            12.6838      235   3497        3732      235
  foldback-lqr-pre-gap                  4229         foldback        0.0051    2228            12.6838      235   497         732       235
  foldback-lqr-start                    4320         foldback        0.0051    2316            11.5801      323   0           323       323
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.3873      0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.3123      0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.3244      0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.3655      0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.3919      0                           0
  large-end-gap-foldback                1720         foldback        0.0119    796             19.3850      0                           0
  longboi                               120000       normal          0.0067                    19.3571      0                           0
  low-quality-first-bases               860          normal          0.0242                    20.5177      0                           0
  low-quality-regions                   712          failed          0.0168                    9.5015       102   443         545       122
  monotonous-read                       60007        failed          0.7926                    19.5242      0                           0
  non-foldback                          2800         normal          0.0073                    19.3556      0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.3504      0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.3448      0                           0
  small-foldback-160bp-total            1160         normal          0.0178                    19.3520      0                           0

  $ minimera --min-qscore 15
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.3244      0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.3460      0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.3582      0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.3441      0                           0
  big-deletion                          3900         normal          0.0052                    19.3885      0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.3682      0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.3935      0                           0
  foldback-longboi                      127949       foldback        0.0072    63970           19.3547      0                           0
  foldback-lqr-end                      4320         failed          0.0051                    11.5801      323   3997        4320      323
  foldback-lqr-mid-gap                  4620         failed          0.0051                    9.3342       626   1997        2623      626
  foldback-lqr-post-gap                 4229         failed          0.0051                    12.6838      235   3497        3732      235
  foldback-lqr-pre-gap                  4229         failed          0.0051                    12.6838      235   497         732       235
  foldback-lqr-start                    4320         failed          0.0051                    11.5801      323   0           323       323
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.3873      0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.3123      0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.3244      0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.3655      0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.3919      0                           0
  large-end-gap-foldback                1720         foldback        0.0119    796             19.3850      0                           0
  longboi                               120000       normal          0.0067                    19.3571      0                           0
  low-quality-first-bases               860          normal          0.0242                    20.5177      0                           0
  low-quality-regions                   712          failed          0.0168                    9.5015       102   443         545       122
  monotonous-read                       60007        failed          0.7926                    19.5242      0                           0
  non-foldback                          2800         normal          0.0073                    19.3556      0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.3504      0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.3448      0                           0
  small-foldback-160bp-total            1160         normal          0.0178                    19.3520      0                           0

Simple vs Dorado qscore computation:

  $ minimera --dorado-mean-qscore
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.3197      0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.3430      0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.3617      0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.3496      0                           0
  big-deletion                          3900         normal          0.0052                    19.3888      0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.3627      0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.4136      0                           0
  foldback-longboi                      127949       foldback        0.0072    63970           19.3546      0                           0
  foldback-lqr-end                      4320         foldback        0.0051    1996            11.5294      323   3997        4320      323
  foldback-lqr-mid-gap                  4620         foldback        0.0051    1993            9.2830       626   1997        2623      626
  foldback-lqr-post-gap                 4229         foldback        0.0051    1993            12.6349      235   3497        3732      235
  foldback-lqr-pre-gap                  4229         foldback        0.0051    2228            12.6349      235   497         732       235
  foldback-lqr-start                    4320         foldback        0.0051    2316            12.2700      323   0           323       323
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.3876      0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.3132      0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.3260      0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.3670      0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.3858      0                           0
  large-end-gap-foldback                1720         foldback        0.0119    796             19.3788      0                           0
  longboi                               120000       normal          0.0067                    19.3572      0                           0
  low-quality-first-bases               860          normal          0.0242                    40.0000      0                           0
  low-quality-regions                   712          normal          0.0168                    9.1560       102   443         545       122
  monotonous-read                       60007        failed          0.7926                    19.5242      0                           0
  non-foldback                          2800         normal          0.0073                    19.3549      0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.3552      0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.3592      0                           0
  small-foldback-160bp-total            1160         normal          0.0178                    19.3502      0                           0

Foldback position epsilon, for end gaps:

  $ minimera --foldback-position-epsilon 140
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.3244      0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.3460      0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.3582      0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.3441      0                           0
  big-deletion                          3900         normal          0.0052                    19.3885      0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.3682      0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.3935      0                           0
  foldback-longboi                      127949       foldback        0.0072    63970           19.3547      0                           0
  foldback-lqr-end                      4320         foldback        0.0051    1996            11.5801      323   3997        4320      323
  foldback-lqr-mid-gap                  4620         foldback        0.0051    1993            9.3342       626   1997        2623      626
  foldback-lqr-post-gap                 4229         foldback        0.0051    1993            12.6838      235   3497        3732      235
  foldback-lqr-pre-gap                  4229         foldback        0.0051    2228            12.6838      235   497         732       235
  foldback-lqr-start                    4320         foldback        0.0051    2316            11.5801      323   0           323       323
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.3873      0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.3123      0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.3244      0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.3655      0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.3919      0                           0
  large-end-gap-foldback                1720         foldback        0.0119    796             19.3850      0                           0
  longboi                               120000       normal          0.0067                    19.3571      0                           0
  low-quality-first-bases               860          normal          0.0242                    20.5177      0                           0
  low-quality-regions                   712          normal          0.0168                    9.5015       102   443         545       122
  monotonous-read                       60007        failed          0.7926                    19.5242      0                           0
  non-foldback                          2800         normal          0.0073                    19.3556      0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.3504      0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.3448      0                           0
  small-foldback-160bp-total            1160         normal          0.0178                    19.3520      0                           0

Large deletion, should be catchable by relaxing intercept epsilon:

  $ minimera --intercept-epsilon 80
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.3244      0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.3460      0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.3582      0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.3441      0                           0
  big-deletion                          3900         foldback        0.0052    1946            19.3885      0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.3682      0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.3935      0                           0
  foldback-longboi                      127949       foldback        0.0072    63970           19.3547      0                           0
  foldback-lqr-end                      4320         foldback        0.0051    1996            11.5801      323   3997        4320      323
  foldback-lqr-mid-gap                  4620         foldback        0.0051    1993            9.3342       626   1997        2623      626
  foldback-lqr-post-gap                 4229         foldback        0.0051    1993            12.6838      235   3497        3732      235
  foldback-lqr-pre-gap                  4229         foldback        0.0051    2228            12.6838      235   497         732       235
  foldback-lqr-start                    4320         foldback        0.0051    2316            11.5801      323   0           323       323
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.3873      0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.3123      0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.3244      0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.3655      0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.3919      0                           0
  large-end-gap-foldback                1720         foldback        0.0119    796             19.3850      0                           0
  longboi                               120000       normal          0.0067                    19.3571      0                           0
  low-quality-first-bases               860          normal          0.0242                    20.5177      0                           0
  low-quality-regions                   712          normal          0.0168                    9.5015       102   443         545       122
  monotonous-read                       60007        failed          0.7926                    19.5242      0                           0
  non-foldback                          2800         normal          0.0073                    19.3556      0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.3504      0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.3448      0                           0
  small-foldback-160bp-total            1160         normal          0.0178                    19.3520      0                           0

Minimum foldback size:

  $ minimera --minimum-cluster-length 20
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.3244      0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.3460      0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.3582      0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.3441      0                           0
  big-deletion                          3900         normal          0.0052                    19.3885      0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.3682      0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.3935      0                           0
  foldback-longboi                      127949       foldback        0.0072    63970           19.3547      0                           0
  foldback-lqr-end                      4320         foldback        0.0051    1996            11.5801      323   3997        4320      323
  foldback-lqr-mid-gap                  4620         foldback        0.0051    1993            9.3342       626   1997        2623      626
  foldback-lqr-post-gap                 4229         foldback        0.0051    1993            12.6838      235   3497        3732      235
  foldback-lqr-pre-gap                  4229         foldback        0.0051    2228            12.6838      235   497         732       235
  foldback-lqr-start                    4320         foldback        0.0051    2316            11.5801      323   0           323       323
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.3873      0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.3123      0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.3244      0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.3655      0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.3919      0                           0
  large-end-gap-foldback                1720         foldback        0.0119    796             19.3850      0                           0
  longboi                               120000       normal          0.0067                    19.3571      0                           0
  low-quality-first-bases               860          normal          0.0242                    20.5177      0                           0
  low-quality-regions                   712          normal          0.0168                    9.5015       102   443         545       122
  monotonous-read                       60007        failed          0.7926                    19.5242      0                           0
  non-foldback                          2800         normal          0.0073                    19.3556      0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.3504      0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.3448      0                           0
  small-foldback-160bp-total            1160         foldback        0.0178    1076            19.3520      0                           0

  $ minimera --minimum-cluster-length 20 --minimum-foldback-length-absolute 200
  $ results
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  481e3d07-9220-434e-9372-d54ff402bf6d  2500         normal          0.0081                    19.3244      0                           0
  6fc48967-2e65-4e92-9c77-6e86e55238d6  2500         normal          0.0079                    19.3460      0                           0
  7e04aba5-e7b1-468a-9319-17110f10f284  2500         normal          0.0099                    19.3582      0                           0
  8d1cadec-dbe5-4668-b9ae-60027aead0a8  2500         normal          0.0080                    19.3441      0                           0
  big-deletion                          3900         normal          0.0052                    19.3885      0                           0
  cae458a3-0575-4f91-870d-8657ba1e74fa  2500         normal          0.0079                    19.3682      0                           0
  cfb7c61a-033b-4085-b79a-eb293b2e1362  500          normal          0.0200                    19.3935      0                           0
  foldback-longboi                      127949       foldback        0.0072    63970           19.3547      0                           0
  foldback-lqr-end                      4320         foldback        0.0051    1996            11.5801      323   3997        4320      323
  foldback-lqr-mid-gap                  4620         foldback        0.0051    1993            9.3342       626   1997        2623      626
  foldback-lqr-post-gap                 4229         foldback        0.0051    1993            12.6838      235   3497        3732      235
  foldback-lqr-pre-gap                  4229         foldback        0.0051    2228            12.6838      235   497         732       235
  foldback-lqr-start                    4320         foldback        0.0051    2316            11.5801      323   0           323       323
  foldback-no-lqr                       4000         foldback        0.0051    1996            19.3873      0                           0
  foldback-tandem-rc-dup                3560         normal          0.0085                    19.3123      0                           0
  full-foldback                         2400         foldback        0.0106    1196            19.3244      0                           0
  gapped-foldback                       2200         foldback        0.0092    1097            19.3655      0                           0
  gapped-partial-foldback               1600         foldback        0.0125    1097            19.3919      0                           0
  large-end-gap-foldback                1720         foldback        0.0119    796             19.3850      0                           0
  longboi                               120000       normal          0.0067                    19.3571      0                           0
  low-quality-first-bases               860          normal          0.0242                    20.5177      0                           0
  low-quality-regions                   712          normal          0.0168                    9.5015       102   443         545       122
  monotonous-read                       60007        failed          0.7926                    19.5242      0                           0
  non-foldback                          2800         normal          0.0073                    19.3556      0                           0
  non-foldback-tandem-rc-dup            2300         normal          0.0108                    19.3504      0                           0
  partial-foldback                      1700         foldback        0.0122    1196            19.3448      0                           0
  small-foldback-160bp-total            1160         normal          0.0178                    19.3520      0                           0

Make sure we can still successfully plot things:

  $ minimera --plot-foldbacks
  $ find out/plots | sort
  out/plots
  out/plots/foldback-longboi.png
  out/plots/foldback-lqr-end.png
  out/plots/foldback-lqr-mid-gap.png
  out/plots/foldback-lqr-post-gap.png
  out/plots/foldback-lqr-pre-gap.png
  out/plots/foldback-lqr-start.png
  out/plots/foldback-no-lqr.png
  out/plots/full-foldback.png
  out/plots/gapped-foldback.png
  out/plots/gapped-partial-foldback.png
  out/plots/large-end-gap-foldback.png
  out/plots/partial-foldback.png

Check that BAM support works:

  $ minimera_on "$TESTDIR/data/ont-giab-hg002.fastq" --fastq --output ./out-fastq
  $ minimera_on "$TESTDIR/data/ont-giab-hg002.bam"   --bam   --output ./out-bam
  $ results_on out-fastq > results_fastq
  $ results_on out-bam > results_bam
  $ cat results_fastq
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  001639d3-c92d-4acd-885c-26d17d48c8e3  9241         normal          0.2201                    17.9220      33    1538        1571      68
  0182bc92-5170-4e2b-8e7d-cb430202f218  7066         normal          0.3124                    13.2011      319   6106        6425      424
  05e47703-9843-4c05-a1dc-01adeec8b597  31124        normal          0.0094                    29.7840      0                           0
  07b3258b-8e7a-4a76-b691-4d6ba1c38835  6653         normal          0.2917                    12.6926      173   5751        5924      294
  08d0c307-91d9-4a07-8879-522094dfa05d  7994         failed          0.6967                    3.6960       1934  2703        4637      5162
  094202bb-1f14-41ed-998c-9d0421168f3d  15174        normal          0.0058                    20.0148      0                           0
  094d7c38-a136-4734-9bca-87171f19f2ce  24765        normal          0.0820                    26.8989      0                           0
  09c63a94-50ed-436f-8ddf-9be58e5a2231  3637         normal          0.3741                    9.7517       299   2781        3080      381
  0ac4669f-775d-42f7-a9ac-80436f1e4829  46159        foldback        0.0157    23071           21.7615      0                           0
  0b4542ad-b8a8-4026-843b-1b82678cfebf  12178        normal          0.0153                    27.7686      0                           0
  0cf1fc94-d6e3-461a-a5ee-fa826be30387  60633        foldback        0.0230    30321           18.1566      0                           0
  0f0a4416-0272-4bb4-81db-e4d258348511  11664        normal          0.0064                    18.0091      0                           0
  109f9b03-32b1-40bd-9407-8cc5d198d25e  3774         normal          0.4651                    12.5908      260   2990        3250      275
  116b375c-5fbc-4410-8044-5076ac45681b  59486        foldback        0.0156    29734           24.4032      0                           0
  12779a77-8412-40cd-8ab3-e1af0cb540fb  34364        normal          0.0115                    27.3530      0                           0
  13ecff2b-2587-4ece-9cda-eaa114bc2fbb  6921         normal          0.0042                    22.1266      0                           0
  16114366-a758-4ea8-9168-5b26333172e3  7013         normal          0.2792                    11.4891      336   5933        6269      526
  16f47723-5bee-4572-a6a6-6d3230406a4f  4514         foldback        0.0140    290             28.5822      0                           0
  17212d5a-08c1-4907-8afe-16c1643c4cfc  17534        normal          0.0164                    23.8327      0                           0
  19ad2153-dbf9-40d7-82cc-d4b695cb3220  6927         foldback        0.0153    3453            22.7259      0                           0
  1b3f2726-36a6-452a-867d-8ef312445350  2628         foldback        0.2966    1321            23.5398      0                           0
  1b936e6e-fa8b-400b-abf9-ee14e2d1ccda  7393         normal          0.3696                    15.2345      243   6384        6627      267
  1bc3d46c-0311-4738-b970-00c0b9e84194  1908         foldback        0.0136    1003            15.9465      0                           0
  1f24e767-1287-4008-af06-c7e399d46b6f  3287         foldback        0.0619    1636            21.6591      0                           0
  21ce0a2c-ae92-48d3-be6d-ff6f6e20042b  448          failed          0.8333                    2.5891       342   106         448       428
  241f5daa-c95f-4a08-b6a2-1a057eeef30d  4947         normal          0.0138                    26.0137      0                           0
  253102b0-f5ea-4e9d-9bcc-21f4d452c745  27174        normal          0.2821                    39.2032      0                           0
  2574d68c-711f-49e1-955d-607a9603726e  6494         normal          0.0574                    27.6294      0                           0
  2678562c-e00f-424a-9cd9-7ef655806ddc  3698         normal          0.4810                    11.3305      295   2873        3168      332
  272ddd5e-976f-465e-b3cc-d3c3653bd76c  11246        foldback        0.0118    5649            16.9089      0                           0
  29277371-b0e8-485b-bbcf-d8bc45ddc63f  44482        foldback        0.0093    22227           20.1518      11    42103       42114     11
  2a7aae42-0b50-46a0-9c0b-7693b69f712a  11385        normal          0.0211                    24.3880      0                           0
  2bce7b82-9969-4eab-b5c5-53f6b0d4585d  47452        foldback        0.0149    23721           19.7415      0                           0
  2c817ff3-6814-47e0-875c-bd11791f2144  43129        foldback        0.0106    21548           15.1783      15    36811       36826     15
  2db84a7e-c9d4-420f-a3ed-e7b9e098be01  29724        normal          0.2708                    26.8544      0                           0
  2e8bbd19-cfe5-4c67-8c32-a529fb6e7866  844          foldback        0.0122    410             17.1992      0                           0
  2fb04cb1-ed1f-4fad-9436-0021bcfc770f  2310         normal          0.0106                    17.7732      0                           0
  2fd370b1-6f03-4423-ad0c-012d59c30b69  333          normal          0.0278                    9.3214       0                           0
  3516fdc8-38f3-4ab9-9025-45ec2555f007  15835        foldback        0.0101    10465           18.0985      0                           0
  354e74fe-b8c6-457e-b7ad-5f8359f37b09  12596        normal          0.0041                    35.1421      0                           0
  37d6d384-4158-4f26-97c4-a61a6515bfc7  15500        normal          0.0061                    24.7484      0                           0
  3a154196-4879-410c-b708-096bff0cb023  607          foldback        0.0328    293             20.9287      0                           0
  3acef803-2277-4bc3-8d5d-418f99673d7f  13756        normal          0.0668                    23.2385      18    0           18        18
  3b90295a-d536-4077-a832-bd01b9a18735  14906        foldback        0.0081    7457            17.5871      0                           0
  3bbb13fc-8257-4fe0-8fe9-004036b5d97d  36035        foldback        0.0139    18069           25.4711      0                           0
  4395dafa-30d4-45b9-bd75-cc2855fb5e4e  3308         normal          0.0154                    29.7423      0                           0
  4401837f-9bc7-4897-9d4d-b17a3deaf979  11887        normal          0.4507                    38.6487      0                           0
  45e468b1-41e5-47ee-befc-3cbe706f071f  6360         foldback        0.0079    3187            14.5524      12    0           12        12
  46e2af21-c709-40f4-812d-11764044fed6  15961        foldback        0.1901    7973            22.5722      0                           0
  47633321-5264-4cfb-9bab-25da442a56f7  20126        normal          0.0100                    22.0998      0                           0
  48a874f0-3e3b-4d88-8880-20990006be84  41598        foldback        0.0123    20792           21.8888      0                           0
  49bfbfe2-cd96-4de1-a322-4ffb08bb5c6b  10256        foldback        0.0039    5120            18.5534      0                           0
  520040a8-525a-4ef4-826c-fa75568afd8b  23597        foldback        0.0159    11809           10.3812      0                           0
  523649be-a9cc-4c58-8cdf-cd0312f301fc  33941        normal          0.0145                    28.8552      0                           0
  52750b61-a1ff-46d9-ac72-eca7321c8d6f  5604         failed          0.2926                    8.8011       1141  4034        5175      1164
  541289bc-75f7-4e24-a498-1a3ed759eaec  15465        foldback        0.1527    7729            20.5239      0                           0
  563d4b0e-9d50-4792-b963-b39468841aa3  18786        foldback        0.0370    9383            21.2787      0                           0
  56f5be8f-091f-4ffe-b942-a9edc6cdab84  9125         normal          0.0060                    13.7235      0                           0
  5cdab79e-fb60-456e-afdf-9a1c04c4a3b3  31988        normal          0.0102                    13.4474      51    14995       15046     87
  5e26595c-6834-4a3d-8026-12b0de365717  30510        normal          0.0113                    29.6000      0                           0
  5e5b8e5a-6b70-4d02-87dd-7bb7c11b9f0e  16023        foldback        0.0110    8005            19.4506      0                           0
  614a6d9f-85d8-4aec-9bcc-9d38f510e419  20557        normal          0.0138                    30.2097      0                           0
  619876b9-278a-4c44-9f58-2b4b114c771e  8716         normal          0.0364                    21.0464      0                           0
  62edc399-b478-41f9-8795-26f44ecf674a  8639         normal          0.2641                    13.5197      331   7762        8093      376
  63caa36b-20c1-4a67-a9f1-64b5e140a82c  6787         normal          0.3521                    13.1220      348   5714        6062      435
  646a4bca-86c5-45b2-b325-018be4a15326  18275        normal          0.0415                    20.9516      0                           0
  64e3094c-388b-416c-bfc0-cbfd3ed9c5ae  23518        foldback        0.0191    11826           18.7727      114   10924       11038     153
  69365b63-3d1d-4ff8-a30a-69692d245cea  7192         normal          0.0166                    21.7972      0                           0
  6a951be4-91d7-4dad-afe3-1d81191736d7  26068        normal          0.2334                    31.1864      0                           0
  73c320df-e7b3-4e7d-896a-4f6574c6bc6d  21931        normal          0.0142                    26.9988      0                           0
  73e14af2-d208-42a8-9f16-24c569f5dd4b  72186        foldback        0.0383    27485           16.6291      159   43532       43691     278
  75d5528f-2764-49a7-8f51-f5e113910495  28441        foldback        0.0529    4509            11.1910      10    18848       18858     10
  77242e55-f94c-435b-9160-9ab8ed916b79  63214        foldback        0.0618    31601           20.8650      0                           0
  77f10c34-4bf4-4d86-a709-bd8a00736942  6033         foldback        0.0066    3028            16.1880      0                           0
  796f0c06-b649-41d6-a33c-ba62bf117bd2  2909         normal          0.0102                    21.3009      0                           0
  7d7af697-7ab6-49cd-9f6c-b33feb1aff15  19842        foldback        0.0135    9913            23.9168      0                           0
  7edc3643-2a77-4b39-bf42-b13aa519bcee  12539        foldback        0.0048    6283            19.1319      0                           0
  7fabf2f2-7fb6-4108-bc3a-24e9aeae1304  32371        normal          0.0953                    34.9136      0                           0
  812215de-bdba-4a4c-9c0a-dcb5ad4a409c  4693         foldback        0.0084    2338            25.4505      0                           0
  853f36b3-60d3-4b5a-b9bc-cb73175e3268  2851         foldback        0.0103    1421            18.6631      0                           0
  894c3fd0-9450-48ca-9fd8-5a67158e836e  5339         failed          0.2663                    5.2772       301   4553        4854      614
  89fb42e0-82f7-423d-a55f-c115b8717f8e  9030         normal          0.0098                    23.9302      0                           0
  8bc1a55a-c07a-4b0f-8824-84d61e7168ed  1001         failed          0.4722                    3.6918       166   835         1001      321
  8ea07fb7-2d49-41ec-9dbf-8e82dea40293  5419         failed          0.8738                    3.0015       879   205         1084      4109
  8f4da53e-ee40-49c4-b529-ab34c2e8fffd  21899        foldback        0.0382    4923            26.0604      0                           0
  90370b6f-ec7f-4882-9f3e-fb2f3da73ec0  11944        normal          0.1637                    24.3657      0                           0
  90d90a38-3a34-48d3-b669-372fcd049025  1193         normal          0.0274                    29.5245      0                           0
  90f7ba19-07ef-4ab2-ac27-3e45a76d5f37  358          normal          0.0274                    28.0363      0                           0
  9155fcee-9842-4916-bac1-bb52ea0d3adf  25825        normal          0.0137                    22.5072      0                           0
  93ff7ccb-d267-4809-a75f-fb9f9943b071  28981        normal          0.3733                    37.0403      0                           0
  949125aa-e812-4d00-94ae-9bda44c5068e  39465        normal          0.0122                    25.2564      0                           0
  9528f0c5-b454-4c53-acf1-89deffe84430  72446        foldback        0.0146    38641           24.8195      0                           0
  957fc627-cb9f-4774-9f7f-9cb539d8a67a  19010        foldback        0.0253    9493            23.9455      0                           0
  967e057e-8050-4521-b80e-23bc7b747b74  7724         normal          0.3443                    12.1001      136   6689        6825      276
  976c9580-2b23-4b08-937a-ba7001eeb4c8  25153        normal          0.0169                    30.0199      0                           0
  97c96fe5-856a-47c8-a604-05548debe63c  21216        normal          0.0063                    30.6748      0                           0
  991afd9d-6ec7-4aa9-b308-76f728e9a8f4  12332        normal          0.0777                    29.8261      0                           0
  99719d6a-d941-4cc2-aebd-db805a76d0b6  6722         normal          0.2795                    12.2008      319   5748        6067      498
  9abdb524-034f-4910-8859-4a3869130d5d  2106         failed          0.9257                    3.6720       594   883         1477      1331
  9b8231dd-ad52-461f-9add-cff77fa13087  38484        foldback        0.0240    19289           10.1705      0                           0
  9bb56f2a-25da-4267-80e0-ca5aeb92c472  15494        normal          0.1914                    27.1429      0                           0
  9cace92b-857a-4568-8ea8-c88d8209f888  17919        foldback        0.0362    15117           13.8956      122   15471       15593     697
  9ddea0ba-04f4-4c51-934d-85a6804a59ff  2463         normal          0.3224                    25.5465      0                           0
  9e2bc8f1-39b9-46b9-82dd-8ca78f15b1cc  9577         normal          0.0294                    33.0734      0                           0
  a19c89e1-3e6f-4a81-9c90-46bf31be985d  12104        normal          0.0148                    9.9898       0                           0
  a22113fe-6b9f-45cf-9710-670b3a8a3a89  1188         foldback        0.0256    603             15.0932      0                           0
  a2e4cf26-ace5-417a-bf7c-e89bb68ef4df  35044        foldback        0.0089    17517           19.8673      0                           0
  a3e0c226-c62b-4fc8-9e03-11aee3fb0971  12811        normal          0.0155                    24.1730      0                           0
  a54279ae-5aad-4efc-9a86-f3b5a7524f32  1274         foldback        0.0154    638             9.8130       0                           0
  a5cb6449-3c51-4c5c-a5c3-1fcb43e928ea  12763        normal          0.0331                    27.9166      0                           0
  a5d4db4d-01af-4428-8d1b-381e3a28d0cf  1268         failed          0.0327                    5.6542       541   0           541       541
  a7d42b5d-1904-4f90-b127-be9e5f9e123d  847          foldback        0.0476    434             16.4957      0                           0
  aa577707-90db-4216-a48c-29a3c6c5ed7d  4163         foldback        0.0187    2073            23.4256      0                           0
  ab378c3e-fe33-41bc-a72e-c950b971bac5  9089         foldback        0.0121    4537            19.5900      0                           0
  ab983b34-0bb9-4caa-b3c7-019b43fbc090  586          failed          0.9914                    6.0208       0                           0
  ac4710bf-6919-479a-aa18-73596a7810ed  31599        normal          0.0146                    26.4584      0                           0
  b1626212-8663-4bf8-a29b-9dc43f75688f  8281         normal          0.0154                    30.3274      0                           0
  b2965558-4b34-49b3-9f1f-916c2da99754  37502        foldback        0.0119    18744           26.7721      0                           0
  b8857dbc-306f-4f73-8354-1c761224829b  13951        normal          0.0382                    29.9175      0                           0
  ba566086-cf90-4814-b8bb-9132f11167c1  18864        normal          0.0130                    21.1147      0                           0
  bc5f2998-25c0-4a13-8fa3-848dc79a17e7  643          foldback        0.0391    315             21.8286      0                           0
  bd8d2e91-1c1f-498c-aeae-74fa27a27180  7076         normal          0.0138                    26.8965      0                           0
  bda9cf3c-a3d4-44c9-9dde-e0acc448070d  5587         normal          0.4185                    13.3573      275   4648        4923      299
  c32da220-3fe5-46fb-860d-9e9bb95bf561  4119         foldback        0.0075    2067            17.4171      0                           0
  c530fc74-10c2-463d-a354-cae070d577be  5568         foldback        0.1736    5164            32.7545      0                           0
  c720befb-f812-4da5-be0a-afc44b6d2727  17446        foldback        0.0100    8724            19.1495      0                           0
  c7ccff9a-bf36-4ac3-9df0-a0a85600ffce  165871       foldback        0.0386    16769           11.1969      763   45002       45765     2555
  c7e3bb5e-c343-4eb3-afb6-186262bd533f  32780        normal          0.4496                    30.9414      0                           0
  c853813d-73d5-4102-8e48-2f300af4dc37  616          foldback        0.0821    402             18.3343      0                           0
  cc6f3714-c02f-432b-9e5b-a7bc4f787d8b  307          normal          0.3333                    41.5691      0                           0
  cd203527-d907-42ea-804d-d969dd219097  11229        foldback        0.0223    10524           31.7994      0                           0
  cf1bf60b-5c22-4a2c-8df1-ecdf60e05b20  13664        normal          0.0188                    27.2365      0                           0
  d044de98-cd4a-4dce-9969-cb0bc5d4bbc4  10411        normal          0.1537                    29.0377      0                           0
  d1226153-a87a-4908-8966-8c56a6db1da4  10626        foldback        0.0079    5324            15.4475      0                           0
  d29bf53a-f739-44ab-bbd4-77ad2428eb5b  5402         normal          0.0674                    25.9278      0                           0
  d4046320-fdec-4d40-8c1d-65c7472dacfc  25196        normal          0.2942                    35.1122      0                           0
  d58f74a1-7794-467e-9d29-ac49d84f0f4b  43           failed          0.5000                    5.9158       10    0           10        10
  d6e815de-ed1f-4a90-9cd5-96ec3b627b94  8331         foldback        0.0055    4160            18.2764      0                           0
  d9a07f01-028c-4e05-8f6d-052e04cc3e3e  31821        normal          0.0135                    26.0763      0                           0
  da16c11e-c9d0-4bfa-8cc8-3a245e8d351f  8524         foldback        0.0075    4255            20.1692      0                           0
  da9606e3-173b-4014-bdff-c85fff5ecfc9  37286        normal          0.1549                    24.6351      0                           0
  de4b2049-60b9-4ef9-a4b1-ee7924fdaffb  12504        foldback        0.0140    398             22.3328      0                           0
  df8b5a14-ade1-4025-8793-513078171550  3679         failed          0.5377                    10.6002      214   2673        2887      287
  df994013-7216-43d1-8101-939caabbcf9e  4113         normal          0.4870                    10.2209      264   3070        3334      368
  e20461e4-f905-40f1-8f45-3c732992070d  22640        normal          0.2439                    33.4333      0                           0
  e2ac6ada-180a-44b9-94c2-a02cfdfe606d  35434        normal          0.0191                    29.5341      0                           0
  e42186d7-1404-474d-8681-b6832e8dc3b4  5400         foldback        0.0065    2717            20.0456      0                           0
  e526dd99-bfd2-4d5a-87bb-55c7879d2666  28011        normal          0.0201                    26.7622      0                           0
  e5f369b5-b19c-44fc-b1f5-c33818c58531  18995        normal          0.0923                    24.7021      0                           0
  e645db33-ac40-4430-b8c7-323b416f275e  28583        normal          0.0244                    29.7844      0                           0
  e7678827-d611-40b2-86be-6a0526a3faee  16339        foldback        0.0185    8162            21.6592      0                           0
  e97ab31e-408b-43ab-8797-456ed449ac06  2047         foldback        0.0277    1798            20.2705      0                           0
  e9e62592-918a-4e27-84b6-f5fe1279610d  7759         normal          0.4618                    12.4995      128   6084        6212      276
  ea9bf365-01bd-4c36-abf8-ae5def9305a7  5961         normal          0.0057                    28.3238      0                           0
  ebe02e28-ad3e-41ad-8ed0-78479942631e  474          foldback        0.0230    228             28.4367      0                           0
  ec528a8d-e692-4d9d-8605-ecedb08613ac  16399        normal          0.1968                    32.0119      0                           0
  ed0f6f04-488e-4825-b535-c7219ada24e0  51995        foldback        0.0099    26007           21.6812      0                           0
  f46f718b-2afb-42ed-874a-b6821488971b  37302        normal          0.0173                    23.7105      0                           0
  f4acee0d-ac37-44ea-b2ef-13992e4cdf7d  82980        foldback        0.0246    41493           21.1314      0                           0
  f9a66946-83d5-4501-b1a5-90f5169c1583  11331        normal          0.0040                    31.6131      0                           0
  fbba2ab0-fc72-4d62-8dde-ce91b271e487  8541         foldback        0.0245    4263            25.4228      0                           0
  fc443da3-857f-4717-a11f-c0671c3c6406  3705         normal          0.3232                    10.4649      453   2872        3325      465
  fce25ed0-0e1b-45d2-aa26-8882b66cbaea  26814        foldback        0.0199    25660           26.9771      0                           0
  fcfabc92-2646-4e54-86ca-73693075c5f4  80958        foldback        0.0305    40501           20.0714      0                           0
  fcfd0c05-1f6e-4523-a968-744e710a985b  6782         normal          0.2779                    13.7047      360   5963        6323      413
  fea68ef0-227f-497b-92e9-e9bf5810232c  86008        foldback        0.0308    43015           20.7004      0                           0
  ff48e816-dd58-4df3-9339-bb87490b0b3f  18954        foldback        0.0159    730             26.8934      0                           0
  $ cat results_bam
  read-id                               read-length  classification  monotony  foldback-point  mean-qscore  llqr  llqr-start  llqr-end  lq-total
  001639d3-c92d-4acd-885c-26d17d48c8e3  9241         normal          0.2201                    17.9220      33    1538        1571      68
  0182bc92-5170-4e2b-8e7d-cb430202f218  7066         normal          0.3124                    13.2011      319   6106        6425      424
  05e47703-9843-4c05-a1dc-01adeec8b597  31124        normal          0.0094                    29.7840      0                           0
  07b3258b-8e7a-4a76-b691-4d6ba1c38835  6653         normal          0.2917                    12.6926      173   5751        5924      294
  08d0c307-91d9-4a07-8879-522094dfa05d  7994         failed          0.6967                    3.6960       1934  2703        4637      5162
  094202bb-1f14-41ed-998c-9d0421168f3d  15174        normal          0.0058                    20.0148      0                           0
  094d7c38-a136-4734-9bca-87171f19f2ce  24765        normal          0.0820                    26.8989      0                           0
  09c63a94-50ed-436f-8ddf-9be58e5a2231  3637         normal          0.3741                    9.7517       299   2781        3080      381
  0ac4669f-775d-42f7-a9ac-80436f1e4829  46159        foldback        0.0157    23071           21.7615      0                           0
  0b4542ad-b8a8-4026-843b-1b82678cfebf  12178        normal          0.0153                    27.7686      0                           0
  0cf1fc94-d6e3-461a-a5ee-fa826be30387  60633        foldback        0.0230    30321           18.1566      0                           0
  0f0a4416-0272-4bb4-81db-e4d258348511  11664        normal          0.0064                    18.0091      0                           0
  109f9b03-32b1-40bd-9407-8cc5d198d25e  3774         normal          0.4651                    12.5908      260   2990        3250      275
  116b375c-5fbc-4410-8044-5076ac45681b  59486        foldback        0.0156    29734           24.4032      0                           0
  12779a77-8412-40cd-8ab3-e1af0cb540fb  34364        normal          0.0115                    27.3530      0                           0
  13ecff2b-2587-4ece-9cda-eaa114bc2fbb  6921         normal          0.0042                    22.1266      0                           0
  16114366-a758-4ea8-9168-5b26333172e3  7013         normal          0.2792                    11.4891      336   5933        6269      526
  16f47723-5bee-4572-a6a6-6d3230406a4f  4514         foldback        0.0140    290             28.5822      0                           0
  17212d5a-08c1-4907-8afe-16c1643c4cfc  17534        normal          0.0164                    23.8327      0                           0
  19ad2153-dbf9-40d7-82cc-d4b695cb3220  6927         foldback        0.0153    3453            22.7259      0                           0
  1b3f2726-36a6-452a-867d-8ef312445350  2628         foldback        0.2966    1321            23.5398      0                           0
  1b936e6e-fa8b-400b-abf9-ee14e2d1ccda  7393         normal          0.3696                    15.2345      243   6384        6627      267
  1bc3d46c-0311-4738-b970-00c0b9e84194  1908         foldback        0.0136    1003            15.9465      0                           0
  1f24e767-1287-4008-af06-c7e399d46b6f  3287         foldback        0.0619    1636            21.6591      0                           0
  21ce0a2c-ae92-48d3-be6d-ff6f6e20042b  448          failed          0.8333                    2.5891       342   106         448       428
  241f5daa-c95f-4a08-b6a2-1a057eeef30d  4947         normal          0.0138                    26.0137      0                           0
  253102b0-f5ea-4e9d-9bcc-21f4d452c745  27174        normal          0.2821                    39.2032      0                           0
  2574d68c-711f-49e1-955d-607a9603726e  6494         normal          0.0574                    27.6294      0                           0
  2678562c-e00f-424a-9cd9-7ef655806ddc  3698         normal          0.4810                    11.3305      295   2873        3168      332
  272ddd5e-976f-465e-b3cc-d3c3653bd76c  11246        foldback        0.0118    5649            16.9089      0                           0
  29277371-b0e8-485b-bbcf-d8bc45ddc63f  44482        foldback        0.0093    22227           20.1518      11    42103       42114     11
  2a7aae42-0b50-46a0-9c0b-7693b69f712a  11385        normal          0.0211                    24.3880      0                           0
  2bce7b82-9969-4eab-b5c5-53f6b0d4585d  47452        foldback        0.0149    23721           19.7415      0                           0
  2c817ff3-6814-47e0-875c-bd11791f2144  43129        foldback        0.0106    21548           15.1783      15    36811       36826     15
  2db84a7e-c9d4-420f-a3ed-e7b9e098be01  29724        normal          0.2708                    26.8544      0                           0
  2e8bbd19-cfe5-4c67-8c32-a529fb6e7866  844          foldback        0.0122    410             17.1992      0                           0
  2fb04cb1-ed1f-4fad-9436-0021bcfc770f  2310         normal          0.0106                    17.7732      0                           0
  2fd370b1-6f03-4423-ad0c-012d59c30b69  333          normal          0.0278                    9.3214       0                           0
  3516fdc8-38f3-4ab9-9025-45ec2555f007  15835        foldback        0.0101    10465           18.0985      0                           0
  354e74fe-b8c6-457e-b7ad-5f8359f37b09  12596        normal          0.0041                    35.1421      0                           0
  37d6d384-4158-4f26-97c4-a61a6515bfc7  15500        normal          0.0061                    24.7484      0                           0
  3a154196-4879-410c-b708-096bff0cb023  607          foldback        0.0328    293             20.9287      0                           0
  3acef803-2277-4bc3-8d5d-418f99673d7f  13756        normal          0.0668                    23.2385      18    0           18        18
  3b90295a-d536-4077-a832-bd01b9a18735  14906        foldback        0.0081    7457            17.5871      0                           0
  3bbb13fc-8257-4fe0-8fe9-004036b5d97d  36035        foldback        0.0139    18069           25.4711      0                           0
  4395dafa-30d4-45b9-bd75-cc2855fb5e4e  3308         normal          0.0154                    29.7423      0                           0
  4401837f-9bc7-4897-9d4d-b17a3deaf979  11887        normal          0.4507                    38.6487      0                           0
  45e468b1-41e5-47ee-befc-3cbe706f071f  6360         foldback        0.0079    3187            14.5524      12    0           12        12
  46e2af21-c709-40f4-812d-11764044fed6  15961        foldback        0.1901    7973            22.5722      0                           0
  47633321-5264-4cfb-9bab-25da442a56f7  20126        normal          0.0100                    22.0998      0                           0
  48a874f0-3e3b-4d88-8880-20990006be84  41598        foldback        0.0123    20792           21.8888      0                           0
  49bfbfe2-cd96-4de1-a322-4ffb08bb5c6b  10256        foldback        0.0039    5120            18.5534      0                           0
  520040a8-525a-4ef4-826c-fa75568afd8b  23597        foldback        0.0159    11809           10.3812      0                           0
  523649be-a9cc-4c58-8cdf-cd0312f301fc  33941        normal          0.0145                    28.8552      0                           0
  52750b61-a1ff-46d9-ac72-eca7321c8d6f  5604         failed          0.2926                    8.8011       1141  4034        5175      1164
  541289bc-75f7-4e24-a498-1a3ed759eaec  15465        foldback        0.1527    7729            20.5239      0                           0
  563d4b0e-9d50-4792-b963-b39468841aa3  18786        foldback        0.0370    9383            21.2787      0                           0
  56f5be8f-091f-4ffe-b942-a9edc6cdab84  9125         normal          0.0060                    13.7235      0                           0
  5cdab79e-fb60-456e-afdf-9a1c04c4a3b3  31988        normal          0.0102                    13.4474      51    14995       15046     87
  5e26595c-6834-4a3d-8026-12b0de365717  30510        normal          0.0113                    29.6000      0                           0
  5e5b8e5a-6b70-4d02-87dd-7bb7c11b9f0e  16023        foldback        0.0110    8005            19.4506      0                           0
  614a6d9f-85d8-4aec-9bcc-9d38f510e419  20557        normal          0.0138                    30.2097      0                           0
  619876b9-278a-4c44-9f58-2b4b114c771e  8716         normal          0.0364                    21.0464      0                           0
  62edc399-b478-41f9-8795-26f44ecf674a  8639         normal          0.2641                    13.5197      331   7762        8093      376
  63caa36b-20c1-4a67-a9f1-64b5e140a82c  6787         normal          0.3521                    13.1220      348   5714        6062      435
  646a4bca-86c5-45b2-b325-018be4a15326  18275        normal          0.0415                    20.9516      0                           0
  64e3094c-388b-416c-bfc0-cbfd3ed9c5ae  23518        foldback        0.0191    11826           18.7727      114   10924       11038     153
  69365b63-3d1d-4ff8-a30a-69692d245cea  7192         normal          0.0166                    21.7972      0                           0
  6a951be4-91d7-4dad-afe3-1d81191736d7  26068        normal          0.2334                    31.1864      0                           0
  73c320df-e7b3-4e7d-896a-4f6574c6bc6d  21931        normal          0.0142                    26.9988      0                           0
  73e14af2-d208-42a8-9f16-24c569f5dd4b  72186        foldback        0.0383    27485           16.6291      159   43532       43691     278
  75d5528f-2764-49a7-8f51-f5e113910495  28441        foldback        0.0529    4509            11.1910      10    18848       18858     10
  77242e55-f94c-435b-9160-9ab8ed916b79  63214        foldback        0.0618    31601           20.8650      0                           0
  77f10c34-4bf4-4d86-a709-bd8a00736942  6033         foldback        0.0066    3028            16.1880      0                           0
  796f0c06-b649-41d6-a33c-ba62bf117bd2  2909         normal          0.0102                    21.3009      0                           0
  7d7af697-7ab6-49cd-9f6c-b33feb1aff15  19842        foldback        0.0135    9913            23.9168      0                           0
  7edc3643-2a77-4b39-bf42-b13aa519bcee  12539        foldback        0.0048    6283            19.1319      0                           0
  7fabf2f2-7fb6-4108-bc3a-24e9aeae1304  32371        normal          0.0953                    34.9136      0                           0
  812215de-bdba-4a4c-9c0a-dcb5ad4a409c  4693         foldback        0.0084    2338            25.4505      0                           0
  853f36b3-60d3-4b5a-b9bc-cb73175e3268  2851         foldback        0.0103    1421            18.6631      0                           0
  894c3fd0-9450-48ca-9fd8-5a67158e836e  5339         failed          0.2663                    5.2772       301   4553        4854      614
  89fb42e0-82f7-423d-a55f-c115b8717f8e  9030         normal          0.0098                    23.9302      0                           0
  8bc1a55a-c07a-4b0f-8824-84d61e7168ed  1001         failed          0.4722                    3.6918       166   835         1001      321
  8ea07fb7-2d49-41ec-9dbf-8e82dea40293  5419         failed          0.8738                    3.0015       879   205         1084      4109
  8f4da53e-ee40-49c4-b529-ab34c2e8fffd  21899        foldback        0.0382    4923            26.0604      0                           0
  90370b6f-ec7f-4882-9f3e-fb2f3da73ec0  11944        normal          0.1637                    24.3657      0                           0
  90d90a38-3a34-48d3-b669-372fcd049025  1193         normal          0.0274                    29.5245      0                           0
  90f7ba19-07ef-4ab2-ac27-3e45a76d5f37  358          normal          0.0274                    28.0363      0                           0
  9155fcee-9842-4916-bac1-bb52ea0d3adf  25825        normal          0.0137                    22.5072      0                           0
  93ff7ccb-d267-4809-a75f-fb9f9943b071  28981        normal          0.3733                    37.0403      0                           0
  949125aa-e812-4d00-94ae-9bda44c5068e  39465        normal          0.0122                    25.2564      0                           0
  9528f0c5-b454-4c53-acf1-89deffe84430  72446        foldback        0.0146    38641           24.8195      0                           0
  957fc627-cb9f-4774-9f7f-9cb539d8a67a  19010        foldback        0.0253    9493            23.9455      0                           0
  967e057e-8050-4521-b80e-23bc7b747b74  7724         normal          0.3443                    12.1001      136   6689        6825      276
  976c9580-2b23-4b08-937a-ba7001eeb4c8  25153        normal          0.0169                    30.0199      0                           0
  97c96fe5-856a-47c8-a604-05548debe63c  21216        normal          0.0063                    30.6748      0                           0
  991afd9d-6ec7-4aa9-b308-76f728e9a8f4  12332        normal          0.0777                    29.8261      0                           0
  99719d6a-d941-4cc2-aebd-db805a76d0b6  6722         normal          0.2795                    12.2008      319   5748        6067      498
  9abdb524-034f-4910-8859-4a3869130d5d  2106         failed          0.9257                    3.6720       594   883         1477      1331
  9b8231dd-ad52-461f-9add-cff77fa13087  38484        foldback        0.0240    19289           10.1705      0                           0
  9bb56f2a-25da-4267-80e0-ca5aeb92c472  15494        normal          0.1914                    27.1429      0                           0
  9cace92b-857a-4568-8ea8-c88d8209f888  17919        foldback        0.0362    15117           13.8956      122   15471       15593     697
  9ddea0ba-04f4-4c51-934d-85a6804a59ff  2463         normal          0.3224                    25.5465      0                           0
  9e2bc8f1-39b9-46b9-82dd-8ca78f15b1cc  9577         normal          0.0294                    33.0734      0                           0
  a19c89e1-3e6f-4a81-9c90-46bf31be985d  12104        normal          0.0148                    9.9898       0                           0
  a22113fe-6b9f-45cf-9710-670b3a8a3a89  1188         foldback        0.0256    603             15.0932      0                           0
  a2e4cf26-ace5-417a-bf7c-e89bb68ef4df  35044        foldback        0.0089    17517           19.8673      0                           0
  a3e0c226-c62b-4fc8-9e03-11aee3fb0971  12811        normal          0.0155                    24.1730      0                           0
  a54279ae-5aad-4efc-9a86-f3b5a7524f32  1274         foldback        0.0154    638             9.8130       0                           0
  a5cb6449-3c51-4c5c-a5c3-1fcb43e928ea  12763        normal          0.0331                    27.9166      0                           0
  a5d4db4d-01af-4428-8d1b-381e3a28d0cf  1268         failed          0.0327                    5.6542       541   0           541       541
  a7d42b5d-1904-4f90-b127-be9e5f9e123d  847          foldback        0.0476    434             16.4957      0                           0
  aa577707-90db-4216-a48c-29a3c6c5ed7d  4163         foldback        0.0187    2073            23.4256      0                           0
  ab378c3e-fe33-41bc-a72e-c950b971bac5  9089         foldback        0.0121    4537            19.5900      0                           0
  ab983b34-0bb9-4caa-b3c7-019b43fbc090  586          failed          0.9914                    6.0208       0                           0
  ac4710bf-6919-479a-aa18-73596a7810ed  31599        normal          0.0146                    26.4584      0                           0
  b1626212-8663-4bf8-a29b-9dc43f75688f  8281         normal          0.0154                    30.3274      0                           0
  b2965558-4b34-49b3-9f1f-916c2da99754  37502        foldback        0.0119    18744           26.7721      0                           0
  b8857dbc-306f-4f73-8354-1c761224829b  13951        normal          0.0382                    29.9175      0                           0
  ba566086-cf90-4814-b8bb-9132f11167c1  18864        normal          0.0130                    21.1147      0                           0
  bc5f2998-25c0-4a13-8fa3-848dc79a17e7  643          foldback        0.0391    315             21.8286      0                           0
  bd8d2e91-1c1f-498c-aeae-74fa27a27180  7076         normal          0.0138                    26.8965      0                           0
  bda9cf3c-a3d4-44c9-9dde-e0acc448070d  5587         normal          0.4185                    13.3573      275   4648        4923      299
  c32da220-3fe5-46fb-860d-9e9bb95bf561  4119         foldback        0.0075    2067            17.4171      0                           0
  c530fc74-10c2-463d-a354-cae070d577be  5568         foldback        0.1736    5164            32.7545      0                           0
  c720befb-f812-4da5-be0a-afc44b6d2727  17446        foldback        0.0100    8724            19.1495      0                           0
  c7ccff9a-bf36-4ac3-9df0-a0a85600ffce  165871       foldback        0.0386    16769           11.1969      763   45002       45765     2555
  c7e3bb5e-c343-4eb3-afb6-186262bd533f  32780        normal          0.4496                    30.9414      0                           0
  c853813d-73d5-4102-8e48-2f300af4dc37  616          foldback        0.0821    402             18.3343      0                           0
  cc6f3714-c02f-432b-9e5b-a7bc4f787d8b  307          normal          0.3333                    41.5691      0                           0
  cd203527-d907-42ea-804d-d969dd219097  11229        foldback        0.0223    10524           31.7994      0                           0
  cf1bf60b-5c22-4a2c-8df1-ecdf60e05b20  13664        normal          0.0188                    27.2365      0                           0
  d044de98-cd4a-4dce-9969-cb0bc5d4bbc4  10411        normal          0.1537                    29.0377      0                           0
  d1226153-a87a-4908-8966-8c56a6db1da4  10626        foldback        0.0079    5324            15.4475      0                           0
  d29bf53a-f739-44ab-bbd4-77ad2428eb5b  5402         normal          0.0674                    25.9278      0                           0
  d4046320-fdec-4d40-8c1d-65c7472dacfc  25196        normal          0.2942                    35.1122      0                           0
  d58f74a1-7794-467e-9d29-ac49d84f0f4b  43           failed          0.5000                    5.9158       10    0           10        10
  d6e815de-ed1f-4a90-9cd5-96ec3b627b94  8331         foldback        0.0055    4160            18.2764      0                           0
  d9a07f01-028c-4e05-8f6d-052e04cc3e3e  31821        normal          0.0135                    26.0763      0                           0
  da16c11e-c9d0-4bfa-8cc8-3a245e8d351f  8524         foldback        0.0075    4255            20.1692      0                           0
  da9606e3-173b-4014-bdff-c85fff5ecfc9  37286        normal          0.1549                    24.6351      0                           0
  de4b2049-60b9-4ef9-a4b1-ee7924fdaffb  12504        foldback        0.0140    398             22.3328      0                           0
  df8b5a14-ade1-4025-8793-513078171550  3679         failed          0.5377                    10.6002      214   2673        2887      287
  df994013-7216-43d1-8101-939caabbcf9e  4113         normal          0.4870                    10.2209      264   3070        3334      368
  e20461e4-f905-40f1-8f45-3c732992070d  22640        normal          0.2439                    33.4333      0                           0
  e2ac6ada-180a-44b9-94c2-a02cfdfe606d  35434        normal          0.0191                    29.5341      0                           0
  e42186d7-1404-474d-8681-b6832e8dc3b4  5400         foldback        0.0065    2717            20.0456      0                           0
  e526dd99-bfd2-4d5a-87bb-55c7879d2666  28011        normal          0.0201                    26.7622      0                           0
  e5f369b5-b19c-44fc-b1f5-c33818c58531  18995        normal          0.0923                    24.7021      0                           0
  e645db33-ac40-4430-b8c7-323b416f275e  28583        normal          0.0244                    29.7844      0                           0
  e7678827-d611-40b2-86be-6a0526a3faee  16339        foldback        0.0185    8162            21.6592      0                           0
  e97ab31e-408b-43ab-8797-456ed449ac06  2047         foldback        0.0277    1798            20.2705      0                           0
  e9e62592-918a-4e27-84b6-f5fe1279610d  7759         normal          0.4618                    12.4995      128   6084        6212      276
  ea9bf365-01bd-4c36-abf8-ae5def9305a7  5961         normal          0.0057                    28.3238      0                           0
  ebe02e28-ad3e-41ad-8ed0-78479942631e  474          foldback        0.0230    228             28.4367      0                           0
  ec528a8d-e692-4d9d-8605-ecedb08613ac  16399        normal          0.1968                    32.0119      0                           0
  ed0f6f04-488e-4825-b535-c7219ada24e0  51995        foldback        0.0099    26007           21.6812      0                           0
  f46f718b-2afb-42ed-874a-b6821488971b  37302        normal          0.0173                    23.7105      0                           0
  f4acee0d-ac37-44ea-b2ef-13992e4cdf7d  82980        foldback        0.0246    41493           21.1314      0                           0
  f9a66946-83d5-4501-b1a5-90f5169c1583  11331        normal          0.0040                    31.6131      0                           0
  fbba2ab0-fc72-4d62-8dde-ce91b271e487  8541         foldback        0.0245    4263            25.4228      0                           0
  fc443da3-857f-4717-a11f-c0671c3c6406  3705         normal          0.3232                    10.4649      453   2872        3325      465
  fce25ed0-0e1b-45d2-aa26-8882b66cbaea  26814        foldback        0.0199    25660           26.9771      0                           0
  fcfabc92-2646-4e54-86ca-73693075c5f4  80958        foldback        0.0305    40501           20.0714      0                           0
  fcfd0c05-1f6e-4523-a968-744e710a985b  6782         normal          0.2779                    13.7047      360   5963        6323      413
  fea68ef0-227f-497b-92e9-e9bf5810232c  86008        foldback        0.0308    43015           20.7004      0                           0
  ff48e816-dd58-4df3-9339-bb87490b0b3f  18954        foldback        0.0159    730             26.8934      0                           0
  $ diff -U0 results_fastq results_bam
