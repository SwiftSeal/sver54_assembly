library(ggplot2)
library(ggtree)
library(dplyr)

# import data ------------------------------------------------------------------

# Resistify extracts multiple NB-ARC domains with _1, _2, notation
# Need a column to express this
resistify <- read.table("../results/resistify/results.tsv",
                        sep = "\t",
                        header = TRUE) %>%
  mutate(gene = gsub("\\.t\\d+$", "", Sequence)) %>%
  mutate(Sequence = paste0(Sequence, "_1"))

rnaseq <- read.table("../results/differential_expression.tsv",
                     sep = "\t",
                     header = TRUE) %>%
  mutate(gene = paste0(gene, ".t1_1")) %>%
  filter(gene %in% resistify$Sequence) %>%
  filter(coef == "Infection")
  

# plot tree --------------------------------------------------------------------

nlr_tree <- read.tree("../results/resistify/nbarc_ced4.tree")

ggplot(rnaseq, aes(x = logFC, colour = coef)) +
  geom_density()

facet_plot(plot,
           data = rnaseq,
           mapping = aes(x = logFC, colour = status),
           panel = "LogFC",
           geom = geom_point) + theme_tree2()

