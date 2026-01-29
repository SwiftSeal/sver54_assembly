library(tidyverse)
library(ggtree)
library(ggthemes)
library(extrafont)
library(scattermore)
#font_import(prompt = FALSE, pattern = "arial")
#loadfonts(device = "win")

theme_set(theme_bw(base_size = 8, base_family = "Arial") + theme(panel.grid = element_blank()))

base_directory = "Z://scratch/sver54_assembly/results"

centromere_positions <- read_tsv(
  "results/centromeres/centromere_positions.bed",
  col_names = FALSE
) |>
  select(chromosome = X1, cen_start = X2, cen_end = X3)

cenh3 <- read_tsv(
  file.path(results_path, "cenh3/trash.regions.bed.gz"),
  col_names = FALSE
) |>
  mutate(length = X3 - X2) |>
  left_join(centromere_positions, by = join_by("X1" == "chromosome")) |>
  mutate(cen_status = ifelse(X2 > cen_start & X3 < cen_end, "inside", "outside"))


cenh3 |>
  filter(length > 500, cen_status == "inside") |>
  ggplot(aes(x = length)) +
  geom_histogram()


ggplot(cenh3 |> filter(cen_status == "inside"), aes(x = length, y = X5, colour = X1)) +
  geom_point(alpha = .2)

cenh3_summary <- cenh3 |>
  group_by(X4) |>
  summarise(mean_coverage = mean(X5), mean_length = mean(length))

ggplot(cenh3, aes(x = mean_length, y = mean_coverage)) +
  geom_point()

ggplot(cenh3 |> filter(X4 == "3254_4"), aes(x = length)) +
  geom_histogram()

cg <- read_tsv("../results/centromeres/centromere_window_CG.bed", col_names = FALSE) |>
  select(chromosome = X1, start = X2, end = X3, CG = X4) |>
  group_by(chromosome) |>
  mutate(end = end - min(start),start = start - min(start))

chg <- read_tsv("../results/centromeres/centromere_window_CHG.bed", col_names = FALSE) |>
  select(chromosome = X1, start = X2, end = X3, CHG = X4) |>
  group_by(chromosome) |>
  mutate(end = end - min(start),start = start - min(start))

chh <- read_tsv("../results/centromeres/centromere_window_CHH.bed", col_names = FALSE) |>
  select(chromosome = X1, start = X2, end = X3, CHH = X4) |>
  group_by(chromosome) |>
  mutate(end = end - min(start),start = start - min(start))

ty3 <- read_tsv("../results/centromeres/centromere_window_ty3.bed", col_names = FALSE) |>
  select(chromosome = X1, start = X2, end = X3, Ty3 = X4) |>
  group_by(chromosome) |>
  mutate(end = end - min(start),start = start - min(start))

trash <- read_tsv("../results/centromeres/centromere_window_trash.bed", col_names = FALSE) |>
  select(chromosome = X1, start = X2, end = X3, TRASH = X7) |>
  group_by(chromosome) |>
  mutate(end = end - min(start), start = start - min(start))

cenh3 <- read_tsv(
  "../results/centromeres/centromere_window_cenh3.tab",
  col_names = FALSE,
  comment = "#"
) |>
  group_by(X1) |>
  mutate(
    CENH3 = (X4 - min(X4))/(max(X4) - min(X4))
  ) |>
  select(chromosome = X1, start = X2, end = X3, CENH3) |>
  group_by(chromosome) |>
  mutate(end = end - min(start), start = start - min(start))

merged <- cg |>
  left_join(chg, by = join_by(chromosome, start, end)) |>
  left_join(chh, by = join_by(chromosome, start, end)) |>
  left_join(ty3, by = join_by(chromosome, start, end)) |>
  left_join(trash, by = join_by(chromosome, start, end)) |>
  left_join(cenh3, by = join_by(chromosome, start, end))

plot <- ggplot(merged) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = CENH3)) +
  scale_fill_gradient(low = "white", high = "#BAB0AC") +
  new_scale_fill() +
  geom_rect(aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0, fill = TRASH)) +
  scale_fill_gradient(low = "white", high = "#76B7B2") +
  new_scale_fill() +
  geom_rect(aes(xmin = start, xmax = end, ymin = 1, ymax = 1.25, fill = Ty3)) +
  scale_fill_gradient(low = "white", high = "#59A14F") +
  geom_line(aes(x = start + (end - start), y = CG), colour = "#4E79A7") +
  geom_line(aes(x = start + (end - start), y = CHG), colour = "#F28E2B") +
  geom_line(aes(x = start + (end - start), y = CHH), colour = "#E15759") +
  facet_grid(chromosome ~ ., scales = "free_x") +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    strip.text.y.right = element_text(angle = 0),
    legend.position = "none"
  ) +
  scale_x_continuous(label = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(
    x = "Centromere length",
    y = "Mean methylation (proportion)"
  )

ggsave("plots/centromere_profile.svg", plot, width = 200, height = 100, units = "mm")
```

```{r}
# import te classifications accoridng to tesorter
tesorter <- read_tsv("../results/tesorter/earlgrey.cls.tsv", skip = 1, col_names = FALSE) |>
  mutate(family = str_to_lower(str_extract(X1, "(rnd-\\d_family-\\d+).*", group = 1))) |>
  select(order = X2, superfamily = X3, clade = X4, family)

# get earlgrey annotations and their positions
earlgrey <- read_tsv(
  "../results/solanum_verrucosum.filteredRepeats.gff", col_names = FALSE
) |>
  filter(!str_detect(X3, "Simple_repeat|Low_complexity")) |>
  filter(!str_detect(X1, "scaffold")) |>
  mutate(family = str_to_lower(str_extract(X9, "ID=(.*?);", group = 1))) |>
  dplyr::select(
    chromosome = X1,
    type = X3,
    start = X4,
    end = X5,
    family
  )

# get centromere positions and join them to earlgrey, identifying if TE lies inside or outside


earlgrey <- earlgrey|>
  left_join(centromere_positions, by = "chromosome") |>
  mutate(
    cen_status = ifelse(start > cen_start & end < cen_end, "inside", "outside")
  )

# create a table of these counts for each family
count_table <- earlgrey |>
  group_by(family, cen_status) |>
  count() |>
  pivot_wider(names_from = cen_status, values_from = n) %>%
  replace(is.na(.), 0)

# values for contingency
total_inside <- sum(count_table$inside)
total_outside <- sum(count_table$outside)

# MATHS!
count_table <- count_table |>
  inner_join(tesorter) |>
  rowwise() |>
  mutate(
    fisher = list(fisher.test(matrix(c(inside, outside, total_inside - inside, total_outside - outside), nrow = 2))),
    p_value = fisher$p.value,
    odds = fisher$estimate
  )

count_table |>
  inner_join(tesorter) |>
  filter(p_value < 0.05) |>
  ggplot(aes(x = odds, fill = clade)) +
  geom_histogram() +
  scale_x_log10()

unusual <- count_table |>
  inner_join(tesorter) |>
  filter(p_value < 0.001) |>
  pull(family)
```

```{r}
tree <- read.tree("../results/tesorter/earlgrey.cls.pep.RT.aln.parstree")
tree$tip.label <- str_extract(tree$tip.label,pattern = "(rnd-\\d_family-\\d+)_.*", group = 1)

plot <- ggtree(tree) %<+% count_table

plot +
  geom_tiplab(align = TRUE, size = 1) +
  geom_tippoint(aes(colour = odds))

plot2 <- ggtree(tidytree::tree_subset(tree, node = "rnd-1_family-236", levels_back = 5)) %<+% count_table
plot2 <- plot2 +
  geom_tippoint(aes(colour = clade)) +
  geom_treescale() +
  scale_colour_tableau()

odds <- count_table |>
  select(family, odds) |>
  column_to_rownames("family")

plot2 <- gheatmap(plot2, odds, color = NA, width = .1) +
  scale_fill_viridis_c(option = "magma")
  
ggsave("plots/centromere_tes.svg", plot2, width = 100, height = 100, units = "mm")
```

# TRASH

trash <- read_csv(
  file.path(results_path, "trash/cpc54.assembly.fa_repeats.csv"),
) |>
  filter(width < 10000) |>
  left_join(centromere_positions, by = join_by("seqID" == "chromosome")) |>
  mutate(cen_status = ifelse(start > cen_start & end < cen_end, "inside", "outside"))

ggplot(trash, aes(x = width)) +
  geom_histogram()

ggplot(trash, aes(x = width, fill = cen_status)) +
  geom_histogram(bins = 200) +
  scale_y_log10()

trash_cenh3 <- read_

# Visualise cen07 TRASH and DANTE repeats

dante <- read_tsv(
  file.path(results_path, "dante/cpc54.dante.gff"),
  col_names = FALSE,
  comment = "#"
) |>
  filter(X4 > 35050000, X5 < 35100000, X1 == "chr07")

trash <- read_csv(
  file.path(results_path, "trash/cpc54.assembly.fa_repeats.csv"),
) |>
  filter(start > 35050000, end < 35100000, seqID == "chr07")

centromere_features <- read_tsv(
  file.path(results_path, "centromeres/centromere_methylation_1kbp.tab"),
  col_names = FALSE,
  comment = "#"
) |>
  filter(X2 > 35050000, X3 < 35100000, X1 == "chr07")

ggplot() +
  geom_rect(data = trash, aes(xmin = start, xmax = end, ymin = 0, ymax = 1), fill = "orange", colour = "white") +
  geom_rect(data = dante, aes(xmin = X4, xmax = X5, ymin = 0, ymax = 1)) +
  geom_line(data = centromere_features, aes(x = X2 + (X3 - X2), y = 1 + X4)) +
  theme_void()

