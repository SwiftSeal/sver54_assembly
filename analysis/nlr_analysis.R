library(ggplot2)
library(ggnewscale)
library(tibble)
library(tidyr)
library(ggtree)
library(dplyr)
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(ComplexHeatmap)
library(viridis)

paired_colours = c(
  "CN" = "#a6cee3", 
  "CNL" = "#1f78b4", 
  "N" = "#b2df8a", 
  "NL" = "#33a02c", 
  "RN" = "#fb9a99", 
  "RNL" = "#e31a1c", 
  "TN" = "#fdbf6f", 
  "TNL" = "#ff7f00"
)
nice_colours = c("#2596be", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#bcbd22", "#17becf", "#7f7f7f")

# import data ------------------------------------------------------------------

resistify <- read.table(
  "../results/resistify/final/results.tsv",
  sep = "\t",
  header = TRUE
) %>%
  mutate(gene = str_replace(Sequence, "\\..*$", "")) %>%
  mutate(helixer = str_detect(gene, "solanum_verrucosum"))

nlr_list <- resistify %>%
  filter(Classification != "None") %>%
  pull(gene)

nlr_list_sequence <- resistify %>%
  filter(Classification != "None") %>%
  pull(Sequence)

lcpm <- read.table("../results/tpm.tsv", header = TRUE) %>%
  mutate(tpm = log10(tpm + 1))

methylation <- read.csv("../results/subfeature_methylation.csv")

homologs <- read.table("../results/refplantnlr_diamond.tsv") %>%
  filter(V3 > 90) %>%
  dplyr::select(Sequence = "V1", Homolog = "V2") %>%
  group_by(Sequence) %>%
  summarise(Homologs = paste(Homolog, collapse = "; "))

de_genes <- read.table("../results/differential_expression.tsv", header = TRUE)

go_terms <- read.table(
  "../results/eggnog/eggnog.emapper.annotations",
  skip = 5,
  sep = "\t",
  fill = TRUE
) %>%
  mutate(gene = gsub(".t\\d", "", V1))

# TE overlap stuff

earlgrey_overlaps <- read.table(
  "../results/final_annotation/final_annotation.longest.earlgrey_overlap.bed"
) %>%
  mutate(proportion = V16 / (V3 - V2)) %>%
  dplyr::select(Sequence = "V4", TE = "V13", proportion) %>%
  filter(proportion > 0.9)

edta_overlaps <- read.table(
  "../results/final_annotation/final_annotation.longest.edta_overlap.bed"
) %>%
  mutate(proportion = V19 / (V3 - V2)) %>%
  filter(proportion > 0.9) %>%
  filter(V12 != "repeat_region") %>%
  dplyr::select(Sequence = "V4", TE = "V12", proportion, V18)

resistify %>%
  filter(Classification != "None") %>%
  inner_join(earlgrey_overlaps) %>%
  group_by(TE) %>%
  summarise(count = n())

resistify %>%
  filter(Classification != "None") %>%
  inner_join(edta_overlaps) %>%
  group_by(TE) %>%
  summarise(count = n())

# Of the helitron overlaps, how many were due to structural annotations?
resistify %>%
  filter(Classification != "None") %>%
  inner_join(edta_overlaps) %>%
  mutate(intact = str_detect(V18, "method=structural")) %>%
  group_by(TE, intact) %>%
  summarise(count = n())

# Just need the classifications for now
methylation <- read.csv("../results/subfeature_methylation.csv") %>%
  mutate(gene = gsub("\\.\\d$", "", gene)) %>%
  filter(column_3 == "exon") %>%
  dplyr::select(gene, classification)

# Big table for all the main numbers

clean_edta <- edta_overlaps %>%
  dplyr::select(Sequence, "EDTA_overlap" = proportion, "EDTA_TE" = TE)

clean_earlgrey <- earlgrey_overlaps %>%
  dplyr::select(Sequence, "Earlgrey_TE" = TE, "Earlgrey_overlap" = proportion)

clean_lcpm <- lcpm %>%
  pivot_wider(names_from = condition, values_from = tpm)

big_table <- resistify %>%
  group_by(gene) %>%
  filter(Length == max(Length)) %>%
  filter(Classification != "None") %>%
  left_join(clean_edta) %>%
  left_join(clean_earlgrey) %>%
  left_join(clean_lcpm) %>%
  left_join(homologs) %>%
  left_join(methylation) %>%
  ungroup() %>%
  mutate(mean_expression = dplyr::select(., "0hr_Pinf_infection":"Temperature_stress_4C") %>% rowMeans())

# Plot tree --------------------------------------------------------------------

nbarc_tree <- read.tree("../results/resistify/final/all_nbarc.msa.raxml.bestTree")
nbarc_tree <- ape::root(nbarc_tree, outgroup = "Ced4")

# Remove the domains from the tree labels
nbarc_tree$tip.label = str_remove(nbarc_tree$tip.label, "_\\d$")

tree_plot <- ggtree(nbarc_tree, layout = "circular") %<+% big_table

tree_plot <- gheatmap(
  tree_plot,
  big_table %>%
    dplyr::select(Sequence, Classification) %>%
    column_to_rownames(var = "Sequence"),
  width = 0.1,
  color = NA,
  colnames = FALSE
) +
  scale_fill_manual(values = paired_colours)

tree_plot <- gheatmap(
  tree_plot + new_scale_fill(),
  big_table %>%
    mutate(
      Feature = case_when(
        MADA == "True" ~ "MADA",
        CJID == "True" ~ "CJID"
      )
    ) %>%
    dplyr::select(Sequence, Feature) %>%
    column_to_rownames(var = "Sequence"),
  width = 0.1,
  offset = 0.5,
  color = NA,
  colnames = FALSE
) +
  scale_fill_manual(values = c(
    "CJID" = nice_colours[7],
    "MADA" = nice_colours[8]
  ), na.value = "white")

te_overlap <- big_table %>%
  dplyr::select(Sequence, EDTA_TE, Earlgrey_TE) %>%
  column_to_rownames(var = "Sequence") %>%
  mutate(
    EDTA_TE = case_when(
      EDTA_TE == "CACTA_TIR_transposon" ~ "TIR",
      EDTA_TE == "Mutator_TIR_transposon" ~ "TIR",
      EDTA_TE == "PIF_Harbinger_TIR_transposon" ~ "TIR",
      EDTA_TE == "Copia_LTR_retrotransposon" ~ "LTR",
      EDTA_TE == "helitron" ~ "Helitron"
    ),
    Earlgrey_TE = case_when(
      Earlgrey_TE == "RC/Helitron" ~ "Helitron",
      Earlgrey_TE == "LTR/Copia" ~ "LTR"
   ))

# Need to simplify te labels

tree_plot <- gheatmap(
  tree_plot + new_scale_fill(),
  te_overlap,
  offset = 1,
  width = 0.2,
  color = NA,
  colnames = FALSE
) +
  scale_fill_manual(
    values = c(
      "Helitron" = nice_colours[4],
      "TIR" = nice_colours[5],
      "LTR" = nice_colours[6]
    ),
    na.value = "white"
  )

ggsave("../../pandoc-thesis/figures/nlr_tree.png", tree_plot + geom_treescale(), width = 5.9, height = 5.9, units = "in")

merged <- big_table %>%
  filter(!is.na(mean_expression)) %>%
  column_to_rownames("Sequence")

homologs_indices <- which(!is.na(merged$Homologs))
homologs_list <- merged$Homologs[!is.na(merged$Homologs)]

ha = rowAnnotation(foo = anno_mark(at = homologs_indices, labels = homologs_list))

ht = draw(
  Heatmap(
    dplyr::select(merged, "0hr_Pinf_infection":"Temperature_stress_4C"),
    km = 20,
    col = magma(100),
    right_annotation = ha,
    show_row_names = FALSE,
  )
)
png("../../pandoc-thesis/figures/nlr_clustered_expression.png", width = 10, height = 10, res = 300, units = "in")
ht
dev.off()

r.dend <- row_dend(ht)
rcl.list <- row_order(ht)

lapply(rcl.list, function(x) length(x))

for (i in 1:length(row_order(ht))) {
  if (i == 1) {
    clu <- t(t(row.names(merged[row_order(ht)[[i]],])))
    out <- cbind(clu, paste0("cluster", i))
    colnames(out) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(merged[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste0("cluster", i))
    out <- rbind(out, clu)
  }
}

# rpi-ver1 analysis ------------------------------------------------------------

ver1_kasp <- read.table(
  "../results/kasp_blast.tsv"
)

annotation_gff <- read.table(
  "../results/final_annotation/final_annotation.gff"
)

# what are the boundaries of the wider locus?
ver1_kasp %>%
  filter(V2 == "chr09") %>%
  filter(V3 == 100) %>%
  summarise(min(V9), max(V9), max(V9) - min(V9))

#and of the smaller?
ver1_kasp %>%
  filter(V2 == "chr09") %>%
  filter(V3 == 100) %>%
  filter(V1 %in% c("DMG400017237", "DMG400017146")) %>%
  summarise(min(V9), max(V9), max(V9) - min(V9))

ver1_kasp %>%
  filter(V2 == "chr09") %>%
  filter(V3 == 100) %>%
  filter(V1 %in% c("DMG400011401", "NLR0226")) %>%
  summarise(min(V9), max(V9), max(V9) - min(V9))

# get genes in smaller locus
small_loci_genes <- annotation_gff %>%
  filter(V3 == "gene") %>%
  filter(V1 == "chr09", V4 > 54606140, V4 < 55951434) %>%
  mutate(gene = gsub("ID=", "", V9)) %>%
  pull(gene)

wider_loci_genes <- annotation_gff %>%
  filter(V3 == "gene") %>%
  filter(V1 == "chr09", V4 > 51196801, V4 < 58124547) %>%
  mutate(gene = gsub("ID=", "", V9)) %>%
  pull(gene)

resistify %>%
  filter(gene %in% wider_loci_genes) %>%
  View()

annotation_gff %>%
  filter(V3 == "gene") %>%
  filter(V1 == "chr09", V4 > 54606140, V4 < 55951434) %>%
  mutate(gene = gsub("ID=", "", V9)) %>%
  left_join(de_genes %>% filter(coef == "Infection")) %>%
  ggplot(aes(xmin = V4, xmax = V5, ymin = 0, ymax = 1, fill = status)) +
  geom_rect()

resistify %>%
  filter(gene %in% small_loci_genes) %>%
  View()
  
go_terms %>%
  filter(gene %in% wider_loci_genes) %>%
  left_join(de_genes, relationship = "many-to-many") %>%
  left_join(resistify, relationship = "many-to-many") %>%
  View()
