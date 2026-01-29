library(tidyverse)

trash_cenh3_depth <- read_tsv(
  "results/cenh3/trash.regions.bed.gz",
  col_names = c("chr", "start", "end", "cenh3_depth")
) |>
  mutate(length = end - start)

centromere_positions <- read_tsv(
  "config/centromere_positions.bed",
  col_names = c("chr", "cen_start", "cen_end")
)

trash_cenh3_depth <- left_join(trash_cenh3_depth, centromere_positions, by = "chr") |>
  mutate(cen_status = ifelse(start > cen_start & end < cen_end, "inside", "outside")) |>
  select(!c("cen_start", "cen_end"))

trash_cenh3_depth |>
  filter(cen_status == "inside") |>
  ggplot(aes(x = length, y = cenh3_depth, colour = cen_status)) +
  geom_point(alpha = 1)

trash_cenh3_depth |>
  ggplot(aes(x = length, colour = cen_status)) +
  geom_histogram() +
  scale_y_log10()
