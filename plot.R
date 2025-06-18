library(tidyverse, quietly=TRUE)

args <- commandArgs(trailingOnly=TRUE)
rlen <- as.numeric(args[1])

# data <- read_csv("hits-seq.csv")
# data <- read_csv("hits-mut.csv")
data <- read_csv(
  "hits.csv",
  col_types = list(col_integer(), col_integer(), col_integer(), col_factor())
)

matches <- data |>
  rename(x=l1, y=l2)

ymax <- max(matches |> pull(y))
xmax <- max(matches |> pull(x))
ymax <- max(c(ymax, xmax))

intercepts <- data |>
  select(c) |>
  mutate(h=sqrt(2 * ((ymax - abs(c))**2))) |>
  mutate(y=c) |>
  mutate(x=0)

intercepts |>
  summarize(cmean=mean(c), csd=sd(c))

intercepts

tic <- 2000

xbreaks <- seq(0, rlen, by=tic)
ybreaks <- c(rev(-seq(0, rlen, by=tic)), seq(tic, rlen, by=tic))

ggplot(matches, aes(x=x, y=y)) +
  coord_fixed(xlim=c(0, rlen),
              ylim=c(-rlen, rlen),
              ratio=1) +
  scale_x_continuous(breaks=xbreaks) +
  scale_y_continuous(breaks=ybreaks) +
  geom_point(shape='o', aes(color=cluster)) +
  geom_segment(x=0,    y=0,    xend=0,    yend=rlen, color="#666666") +
  geom_segment(x=rlen, y=0,    xend=rlen, yend=rlen, color="#666666") +
  geom_segment(x=0,    y=0,    xend=rlen, yend=0, color="#666666") +
  geom_segment(x=0,    y=rlen, xend=rlen, yend=rlen, color="#666666") +
  geom_segment(x=0,    y=0,    xend=rlen, yend=rlen, color="#AA6666", linetype="dashed") +
  geom_jitter(color="#00000066", shape='o', data=intercepts) +
  labs(
    title="Mimimizer matches",
    x="Location in Read 1",
    y="Location in Read 2",
  )

ggsave("matches.png", width=6, height=8, unit="in", dpi=600)


ggplot(intercepts, aes(x=y)) +
  # geom_histogram(bins=100, aes(y=after_stat(count / pmax(20000, sqrt(2 * ((ymax - abs(x))**2)))))) +
  geom_histogram(bins=30) +
  labs(
    title="y-intercepts of matches",
  ) + geom_line(aes(y=h/100000000))

ggsave("intercepts.png")
