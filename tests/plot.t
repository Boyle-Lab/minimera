  $ . "$TESTDIR/setup.sh"

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
