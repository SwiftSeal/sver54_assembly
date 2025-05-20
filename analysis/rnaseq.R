
library(DESeq2)
library(tximport)
library(tidyverse)
library(ggpubr)
library(ggthemes)

get_contrast <- function(dds, coef) {
  base_condition <- str_split_i(coef, "_vs_", -1)
  dds$condition <- relevel(dds$condition, base_condition)
  dds <- nbinomWaldTest(dds)
  res <- lfcShrink(dds, coef = coef, type = "apeglm")
  return (res)
}

## Import data

tx2gene <- read_tsv("../results/nf_rnaseq/tx2gene.tsv", col_names = FALSE)

samples <- read_csv("../config/rnaseq.csv") |>
  mutate(condition = str_remove(sample, "_REP\\d")) |>
  mutate(quant = fs::path("../results/nf_rnaseq", sample, "quant.sf"))

txi <- tximport(
  samples$quant,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

ddsTxi <- DESeqDataSetFromTximport(
  txi, 
  colData = samples,
  design = ~ condition
)

dds <- DESeq(ddsTxi)

vsd <- vst(dds)

pca_data <- DESeq2::plotPCA(vsd, returnData = TRUE)

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point() +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank()) +
  scale_color_tableau(
    name = NULL,
    labels = c(
      "INFECTION_0" = "Infection 0hpi",
      "INFECTION_24" = "Infection 24hpi",
      "LEAF" = "Leaf",
      "ROOT" = "Root",
      "SHOOT" = "Shoot",
      "TEMPERATURE_25" = "25°C",
      "TEMPERATURE_35" = "35°C",
      "TEMPERATURE_4" = "4°C"
    )
  )

pca_plot

normalised_counts <- counts(dds, normalized = TRUE) |>
  as.data.frame()

colnames(normalised_counts) <- dds$sample
normalised_counts$gene <- rownames(normalised_counts)

normalised_counts <- normalised_counts |>
  pivot_longer(!gene, names_to = "sample") |>
  mutate(condition = str_remove(sample, "_REP\\d")) |>
  group_by(gene, condition) |>
  summarise(value = mean(value))

# Experimental analysis --------------------------------------------------------
# I've split the RNAseq analysis into individual experiments, as some samples
# can't be compared directly due to tissue differences.

infection_samples <- samples |>
  filter(condition %in% c("INFECTION_24", "INFECTION_0"))

txi <- tximport(
  infection_samples$quant,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

ddsTxi <- DESeqDataSetFromTximport(
  txi, 
  colData = infection_samples,
  design = ~ condition
)

dds <- DESeq(ddsTxi)

res <- get_contrast(dds, "condition_INFECTION_24_vs_INFECTION_0")
summary(res)

infection_table <- as.data.frame(res) |>
  mutate(contrast = "Infection")

tissue_samples <- samples |>
  filter(condition %in% c("LEAF", "SHOOT", "ROOT"))

txi <- tximport(
  tissue_samples$quant,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

ddsTxi <- DESeqDataSetFromTximport(
  txi, 
  colData = tissue_samples,
  design = ~ condition
)

dds <- DESeq(ddsTxi)

res <- get_contrast(dds, "condition_ROOT_vs_LEAF")
root_leaf_table <- as.data.frame(res) |>
  mutate(contrast = "root_vs_leaf")

res <- get_contrast(dds, "condition_SHOOT_vs_LEAF")
shoot_leaf_table <- as.data.frame(res) |>
  mutate(contrast = "shoot_vs_leaf")

res <- get_contrast(dds, "condition_SHOOT_vs_ROOT")
shoot_root_table <- as.data.frame(res) |>
  mutate(contrast = "shoot_vs_root")

temperature_samples <- samples |>
  filter(condition %in% c("TEMPERATURE_25", "TEMPERATURE_4", "TEMPERATURE_35"))

txi <- tximport(
  temperature_samples$quant,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

ddsTxi <- DESeqDataSetFromTximport(
  txi, 
  colData = temperature_samples,
  design = ~ condition
)

dds <- DESeq(ddsTxi)

res <- get_contrast(dds, "condition_TEMPERATURE_4_vs_TEMPERATURE_25")
cold_table <- as.data.frame(res) |>
  mutate(contrast = "cold")

res <- get_contrast(dds, "condition_TEMPERATURE_35_vs_TEMPERATURE_25")
heat_table <- as.data.frame(res) |>
  mutate(contrast = "heat")

de_table <- rbind(
  infection_table,
  root_leaf_table,
  shoot_leaf_table,
  shoot_root_table,
  cold_table,
  heat_table
) |>
  mutate(status = case_when(
    padj < 0.01 & log2FoldChange > 1 ~ "up-regulated",
    padj < 0.01 & log2FoldChange < -1 ~ "down-regulated",
    .default = "unchanged"
  ))

de_plot <- de_table |>
  filter(status != "unchanged") |>
  ggplot(aes(y = contrast, fill = status)) +
  geom_bar(position = "dodge", colour = "#333333") +
  labs(x = "# differentially expressed genes", y = NULL) +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank()) +
  scale_y_discrete(
    labels = c(
      "shoot_vs_root" = "Shoot vs. root",
      "shoot_vs_leaf" = "Shoot vs. leaf",
      "root_vs_leaf" = "Root vs. leaf",
      "Infection" = "Infection",
      "heat" = "Heat stress",
      "cold" = "Cold stress"
    )  
  ) +
  scale_fill_tableau(name = NULL)

final_plot <- ggarrange(
  pca_plot,
  de_plot,
  ncol = 1,
  align = "hv",
  labels = "auto",
  font.label = list(size = 10)
)

write_csv(normalised_counts, file = "data/rnaseq_normalised_counts.csv")
write_csv(de_table, file = "data/rnaseq_de_table.csv")
ggsave("plots/rnaseq.png", final_plot, width = 5.9, height = 4, dpi = 600)