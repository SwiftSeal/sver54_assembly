library(ggtree)
library(dplyr)

tree <- read.tree("../results/tesorter/tesorter.cls.pep.RT.aln.tree")
tip_labels <- tree$tip.label
tip_labels <- sub("#.*", "", tip_labels)
tip_labels <- sub("'", "", tip_labels)
tree$tip.label <- tip_labels

tesorter <- read.table(
  "../results/tesorter/tesorter.cls.tsv", 
  skip = 1,
  comment.char = "",
  sep = "\t"
  ) %>%
  mutate(V1 = sub("#.*", "", V1))
  

p <- ggtree(tree, layout = "equal_angle") %<+% tesorter
p + geom_tippoint(aes(colour = V3))
