library(tidyverse)
library(ggplot2)

data <- read.csv("data/phred.downsampled.csv", col.names=c("x", "qs", "qad", "qmd"))

data <- data |>
  filter(x >= -4000) |>
  filter(x <= 12000) |>
  mutate(xf=as.factor((x %/% 100) * 100))

# Dot plot of quality scores
p <- ggplot(data) +
  geom_jitter(aes(x=x, y=qs), alpha=0.008, width=0, height=0.2, shape='o') +
  labs(title="Plot of Q scores",
       x="Position in read relative to foldback point",
       y="Q score")
ggsave('plot0.png', units='px', plot=p, width=2400, height=1300, dpi=160)

# Bin plot of quality scores
p <- ggplot(data) +
  geom_bin_2d(aes(x=x, y=qs), binwidth=c(100, 1)) +
  scale_fill_gradient(low="white", high="black") +
  labs(title="Bin plot of Q scores",
       x="Position in read relative to foldback point",
       y="Q score")
ggsave('plot1.png', units='px', plot=p, width=2400, height=1300, dpi=160)

# Bin plot of quality - mean scores
p <- ggplot(data) +
  geom_bin_2d(aes(x=x, y=qad), binwidth=c(100, 1)) +
  scale_fill_gradient(low="white", high="black") +
  labs(title="Bin plot of Q - mean Q scores",
       x="Position in read relative to foldback point",
       y="Q score - mean Q score of the read")
ggsave('plot2.png', units='px', plot=p, width=2400, height=1300, dpi=160)

# Bin plot of quality - median scores
p <- ggplot(data) +
  geom_bin_2d(aes(x=x, y=qmd), binwidth=c(100, 1)) +
  scale_fill_gradient(low="white", high="black") +
  labs(title="Bin plot of Q - median Q scores",
       x="Position in read relative to foldback point",
       y="Q score - median Q score of the read")
ggsave('plot3.png', units='px', plot=p, width=2400, height=1300, dpi=160)

# Tukey plots of quality - median scores
p <- ggplot(data) +
  geom_boxplot(aes(x=xf, y=qmd), orientation='x', outlier.alpha=0.005) +
  labs(title="Tukey plots of Q - median Q scores",
       x="Position in read relative to foldback point",
       y="Q score - median Q score of the read")
ggsave('plot4.png', units='px', plot=p, width=2400, height=1300, dpi=160)
