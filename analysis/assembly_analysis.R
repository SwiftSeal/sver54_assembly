library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# import data ------------------------------------------------------------------

gc <- read.table("../results/windows/gc.bed") %>%
  select(V1, V2, V3, V5) %>%
  rename(chrom = "V1", start = "V2", end = "V3", prop = "V5") %>%
  mutate(feature = "GC")

genes <- read.table("../results/windows/genes.bed") %>%
  select(V1, V2, V3, V7) %>%
  rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
  mutate(feature = "Gene")

copia <- read.table("../results/windows/copia.bed") %>%
  select(V1, V2, V3, V7) %>%
  rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
  mutate(feature = "Ty1-copia")

gypsy <- read.table("../results/windows/gypsy.bed") %>%
  select(V1, V2, V3, V7) %>%
  rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
  mutate(feature = "Ty3")

tir <- read.table("../results/windows/tir.bed") %>%
  select(V1, V2, V3, V7) %>%
  rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
  mutate(feature = "TIR")

helitron <- read.table("../results/windows/helitron.bed") %>%
  select(V1, V2, V3, V7) %>%
  rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
  mutate(feature = "Helitron")

methylation <- read.table("../results/windows/methylation.tab") %>%
  rename(
    chrom = "V1",
    start = "V2",
    end = "V3",
    CG = "V4",
    CHG = "V5",
    CHH = "V6"
  ) %>%
  pivot_longer(cols = CG:CHH, names_to = "feature", values_to = "prop")

merged <- rbind(gc, genes, copia, gypsy, tir, helitron, methylation)

plot <- ggplot(
    merged,
    aes(xmin = start, xmax = end, ymin = 0, ymax = prop, fill = prop)
  ) +
  geom_rect() +
  facet_grid(factor(feature, levels = c("GC", "CG", "CHG", "CHH", "Gene", "Ty3", "Ty1-copia", "TIR", "Helitron")) ~ chrom, scales = "free") +
  theme(
    #panel.spacing = unit(0, "lines"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text.y.right = element_text(angle = 0),
    legend.position = "none"
  ) +
  scale_fill_viridis(option = "rocket")
plot

ggsave("genome_features.png", width = 10, height = 6, units = "in", dpi = 600)
