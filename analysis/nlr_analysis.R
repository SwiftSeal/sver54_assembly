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
nice_colours = c("#2596be", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

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

tpm <- read.table("../results/tpm.tsv") %>%
  mutate(tpm = log2(tpm + 1))

tpm_matrix <- tpm %>%
  inner_join(resistify, relationship = "many-to-many") %>%
  dplyr::select(Sequence, condition, tpm) %>%
  pivot_wider(names_from = condition, values_from = tpm) %>%
  column_to_rownames(var = "Sequence")

methylation <- read.table(
  "../results/deepsignal/gene_matrix.tab",
  comment.char = "@"
) %>%
  dplyr::select(!c(V1, V2, V3, V5, V6)) %>%
  filter(V4 %in% nlr_list)

colnames(methylation) <- c(
  "gene",
  paste0("A", 1:600),
  paste0("B", 1:600),
  paste0("C", 1:600)
)

methylation <- methylation %>%
  pivot_longer(!gene, names_to = "origin", values_to = "value") %>%
  separate(origin, into = c("Sample", "Position"), sep = 1)

summary_methylation <- methylation %>%
  group_by(Sample, Position) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    stderr_min = mean(value, na.rm = TRUE) - sd(value, na.rm = TRUE)/sqrt(n()),
    stderr_max = mean(value, na.rm = TRUE) + sd(value, na.rm = TRUE)/sqrt(n()),
  )

ggplot(summary_methylation, aes(x = as.numeric(Position), ymin = stderr_min, ymax = stderr_max, fill = Sample)) +
  geom_line(aes(y = mean)) +
  geom_ribbon(aes(alpha = 0.5))

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
  View()
  summarise(count = n())


# Big table for all the main numbers

clean_edta <- edta_overlaps %>%
  dplyr::select(Sequence, "EDTA_overlap" = proportion, "EDTA_TE" = TE)

clean_earlgrey <- earlgrey_overlaps %>%
  dplyr::select(Sequence, "Earlgrey_TE" = TE, "Earlgrey_overlap" = proportion)

clean_tpm <- tpm %>%
  pivot_wider(names_from = condition, values_from = tpm)

low_expression <- tpm %>%
  group_by(gene) %>%
  summarise(tpm = mean(tpm)) %>%
  mutate("Low expression" = tpm < 2.5) %>%
  dplyr::select(gene, "Low expression")

big_table <- resistify %>%
  group_by(gene) %>%
  filter(Length == max(Length)) %>%
  filter(Classification != "None") %>%
  left_join(clean_edta) %>%
  left_join(clean_earlgrey) %>%
  left_join(clean_tpm) %>%
  left_join(low_expression) %>%
  left_join(homologs) %>%
  ungroup()

ver1_kasp <- read.table("../results/kasp_blast.tsv")

annotation_gff <- read.table("../results/final_annotation/final_annotation.gff")

## Summarise NLRs --------------------------------------------------------------

resistify %>%
  filter(Classification != "None") %>%
  summarise(count = n())

resistify %>%
  group_by(Classification) %>%
  summarise(count = n())

# Plot tree --------------------------------------------------------------------

tree_plot <- ggtree(nbarc_tree) %<+% big_table

classifications <- big_table %>%
  dplyr::select(Sequence, Classification) %>%
  column_to_rownames(var = "Sequence")

tree_plot <- gheatmap(
  tree_plot,
  classifications,
  width = 0.1,
  color = NA) +
  scale_fill_manual(values = paired_colours)

Feature <- big_table %>%
  mutate(Feature = case_when(
    MADA == "True" ~ "MADA",
    CJID == "True" ~ "CJID"
    )
  ) %>%
  dplyr::select(Sequence, Feature) %>%
  column_to_rownames(var = "Sequence")

tree_plot <- gheatmap(
  tree_plot + new_scale_fill(),
  Feature,
  width = 0.1,
  offset = 1,
  color = NA,
  colnames_angle = 90
)

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
  offset = 2,
  width = 0.2,
  color = NA,
  colnames_angle = 90
)

tree_plot <- gheatmap(
  tree_plot + new_scale_fill(),
  tpm_matrix,
  offset = 4,
  color = NA,
  colnames_angle = 90
) +
  scale_fill_viridis_c(option = "magma")

ggsave("nlr_tree.pdf", tree_plot + ylim(-100, 600), width = 8, height = 6, units = "in")

# Expression analysis ----------------------------------------------------------

de_genes %>%
  filter(gene %in% nlr_list) %>%
  filter(status != "not significant") %>%
  ggplot(aes(x = coef, fill = status)) +
  geom_bar()

joined <- inner_join(resistify, tpm, relationship = "many-to-many") %>%
  filter(Classification != "None")

joined %>%
  ggplot(aes(y = tpm, fill = condition)) +
  geom_boxplot()

# Are there differences in NLR expression between conditions?
condition_lm <- summary(lm(tpm ~ condition, data = joined))
condition_lm

expression_histogram <- joined %>%
  ggplot(aes(x = tpm, fill = condition)) +
  geom_histogram() +
  labs(x = "Log2(TPM)", y = "# NLRs") +
  scale_fill_manual(values = nice_colours)
expression_histogram

# how many are low to unexpressed?

joined %>%
  group_by(Sequence) %>%
  summarise(tpm = mean(tpm)) %>%
  summarise(sum(tpm < 1) / n())
# 23.1%

# Is this cus of helixer?

helixer_lm <- summary(lm(tpm ~ helixer, data = joined))

helixer_resistify <- helixer_resistify %>%
  mutate(Sequence = gsub("\\.1", "", Sequence)) %>%
  mutate(Missing = Sequence %in% missing_nlrs$gene_id)


nlr_clustering <- data.frame(Type = character(), Count = integer())

for(i in seq(0, 100000, by = 1000)) {
  overlaps <- countOverlaps(genes(helixer_gff)[helixer_nlrs,], maxgap = i)
  print(sum(overlaps == 1))
  print(sum(overlaps > 1))
}

helixer_resistify <- helixer_resistify %>%
  mutate(Overlapping = Sequence %in% earlgrey_intersect$V4)

tree_plot <- ggtree(helixer_tree, layout = "daylight") %<+% helixer_resistify
tree_plot +
  geom_tippoint(aes(colour = Missing, shape = Overlapping))

# What NLRs overlap with EarlGrey TEs?


merged <- inner_join(resistify, tpm, relationship = "many-to-many") %>%
  filter(Classification != "None") %>%
  dplyr::select(Sequence, condition, tpm) %>%
  pivot_wider(names_from = "condition", values_from = "tpm") %>%
  column_to_rownames("Sequence")

ht = draw(Heatmap(merged, km = 20, col = magma(100)))
pdf("nlr_clustered.pdf", width = 20, height = 20)
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
