  $ . "$TESTDIR/setup.sh"

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
