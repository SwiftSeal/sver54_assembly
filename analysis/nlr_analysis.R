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
                     header = TRUE)

# plot tree --------------------------------------------------------------------

nlr_tree <- read.tree("../results/resistify/nbarc_ced4.tree")

ggtree(nlr_tree) %<+% resistify +
  geom_tippoint(aes(colour = Classification))
