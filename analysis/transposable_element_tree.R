library(ggtree)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(treeio)

tree <- read.tree(
    "results/tesorter/earlgrey.cls.pep.RT.aln.parstree"
)

classifications <- as.data.frame(tree$tip.label) %>%
    mutate(
        id = str_extract(tree$tip.label, "^(rnd-\\d+_family-\\d+)_(.*?)(?:/(.*?))?(?:/(.*?))?$", 1),
        type = str_extract(tree$tip.label, "^(rnd-\\d+_family-\\d+)_(.*?)(?:/(.*?))?(?:/(.*?))?$", 2),
        family = str_extract(tree$tip.label, "^(rnd-\\d+_family-\\d+)_(.*?)(?:/(.*?))?(?:/(.*?))?$", 3),
        subfamily = str_extract(tree$tip.label, "^(rnd-\\d+_family-\\d+)_(.*?)(?:/(.*?))?(?:/(.*?))?$", 4)
    )


tree$tip.label <- classifications$id

centromere_association <- read.csv(
    "results/earlgrey_centromere_enrichment.csv"
) %>%
    mutate(id = tolower(id)) %>%
    filter(id %in% tree$tip.label) %>%
    left_join(classifications, by = "id")

# Might as well summarise here as well eh
centromere_association %>%
    filter(p_value < 0.05 & odds_ratio > 1) %>%
    group_by(type, subfamily) %>%
    summarise(
        n = n(),
        odds_ratio = max(odds_ratio)
    ) %>%
    arrange(odds_ratio, desc(odds_ratio))

centromere_association %>%
    arrange(odds_ratio)

plot_tree <- ggtree(tree, layout = "daylight")

zoomed <- tidytree::tree_subset(tree, "rnd-1_family-627", levels_back = 9)

zoomed_plot <- ggtree(zoomed) %<+% centromere_association +
    geom_tiplab(aes(label = label), size = 2, offset = 0.002) +
    geom_tippoint(aes(colour = odds_ratio)) +
    scale_colour_gradient2(midpoint = 10)

save_plot <- plot_tree +
    geom_tippoint(aes(subset = isTip & label == "rnd-1_family-627"), colour = "red", shape = 21, size = 5) +
    theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
    )


save_plot <- zoomed_plot + patchwork::inset_element(save_plot, left = -0.05, bottom = 0.55, right = 0.35, top = 1.2)

ggsave("../pandoc-thesis/figures/centromere_association_tree.png", save_plot, width = 5.9, height = 5, units = "in", dpi = 600)















