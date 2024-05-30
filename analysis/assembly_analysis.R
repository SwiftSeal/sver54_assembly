library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(Biostrings)
library(viridis)

chromosomes <- c(
  "chr01",
  "chr02",
  "chr03",
  "chr04",
  "chr05",
  "chr06",
  "chr07",
  "chr08",
  "chr09",
  "chr10",
  "chr11",
  "chr12"
)

# import data ------------------------------------------------------------------

gc <- read.table("../results/windows/gc.bed") %>%
  dplyr::select(V1, V2, V3, V5) %>%
  dplyr::rename(chrom = "V1", start = "V2", end = "V3", prop = "V5") %>%
  mutate(feature = "GC")

genes <- read.table("../results/windows/genes.bed") %>%
  dplyr::select(V1, V2, V3, V7) %>%
  dplyr::rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
  mutate(feature = "Gene")

ty1 <- read.table("../results/windows/Ty1.bed") %>%
  dplyr::select(V1, V2, V3, V7) %>%
  dplyr::rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
  mutate(feature = "Ty1")

ty3 <- read.table("../results/windows/Ty3.bed") %>%
  dplyr::select(V1, V2, V3, V7) %>%
  dplyr::rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
  mutate(feature = "Ty3")

cenh3 <- read.csv("../results/windows/cenh3.bed") %>%
  filter(V1 %in% chromosomes) %>%
  mutate(prop = V1.1 / max(V1.1), feature = "CENH3") %>%
  dplyr::select(chrom = "V1", start = "V2", end = "V3", prop, feature)

#tir <- read.table("../results/windows/tir.bed") %>%
#  select(V1, V2, V3, V7) %>%
#  rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
#  mutate(feature = "TIR")
#
#helitron <- read.table("../results/windows/helitron.bed") %>%
#  select(V1, V2, V3, V7) %>%
#  rename(chrom = "V1", start = "V2", end = "V3", prop = "V7") %>%
#  mutate(feature = "Helitron")
#

methylation <- read.table("../results/windows/methylation.tab") %>%
  dplyr::rename(
    chrom = "V1",
    start = "V2",
    end = "V3",
    CG = "V4",
    CHG = "V5",
    CHH = "V6"
  ) %>%
  pivot_longer(cols = CG:CHH, names_to = "feature", values_to = "prop")

merged <- rbind(gc, genes, ty1, ty3, cenh3, methylation)



merged <- merged %>%
  filter(chrom %in% chromosomes)

plot <- ggplot(
    merged,
    aes(xmin = start, xmax = end, ymin = 0, ymax = prop, fill = prop)
  ) +
  geom_rect() +
  facet_grid(factor(feature, levels = c("GC", "CG", "CHG", "CHH", "CENH3", "Gene", "Ty3", "Ty1")) ~ chrom, scales = "free") +
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
  scale_fill_viridis(option = "magma")

ggsave(
  "genome_features.png",
  plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600
)
