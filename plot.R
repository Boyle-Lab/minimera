library(tidyverse)

# data <- read_csv("hits-seq.csv")
data <- read_csv("hits-mut.csv")

matches <- data |>
  rename(x=l1, y=l2)

intercepts <- data |>
  select(c) |>
  rename(y=c) |>
  mutate(x=0)


ggplot(matches, aes(x=x, y=y)) +
  geom_point(shape='o') +
  geom_jitter(color="#FF000066", shape='o', data=intercepts, width=5, height=1) +
  labs(
    title="Mimimizer matches",
    x="Location in Read 1",
    y="Location in Read 2",
  )


# ggplot(intercepts, aes(x=y)) +
#   geom_histogram() +
#   labs(
#     title="y-intercepts of matches",
#   )
