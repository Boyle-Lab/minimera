  $ . "$TESTDIR/setup.sh"

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
