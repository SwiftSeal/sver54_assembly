library(dplyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(limma)
library(edgeR)
library(tximport)

# Import data ------------------------------------------------------------------

# Experiment metadata (samples, repeats, etc)
metadata <- read.table("../config/rna_seq_samples.tsv", header = TRUE)

# Get transcript to gene relationship table from gff file
transcript_to_gene <- read.table(
  "../results/solanum_verrucosum.gtf",
  skip = 2,
  sep = "\t"
  ) %>%
  filter(V3 == "transcript") %>%
  mutate(
    gene = str_extract(V9, "gene_id (\\w+);", group = 1),
    transcript = str_extract(V9, "transcript_id ([^;]+);", group = 1)
  ) %>%
  select(c(transcript, gene))
  
# Salmon quantification files
salmon_quantification <- tximport(
  file.path("../results/salmon/quant", metadata$sample, "quant.sf"),
  type = "salmon",
  tx2gene = transcript_to_gene,
  countsFromAbundance = "lengthScaledTPM"
)

# Analysis ---------------------------------------------------------------------

design_matrix <- model.matrix(~ 0 + condition, data = metadata)

# make contrasts
contrast_matrix <- makeContrasts(
  Infection = conditionInfection_24 - conditionInfection_0,
  RootvsLeaf = conditionRoot - conditionLeaf,
  RootvsShoot = conditionRoot - conditionShoot,
  ShootvsLeaf = conditionShoot - conditionLeaf,
  Temperature35vs25 = conditionTemperature_35C - conditionTemperature_25C,
  Temperature35vs4 = conditionTemperature_35C - conditionTemperature_4C,
  Temperature25vs4 = conditionTemperature_25C - conditionTemperature_4C,
  levels = colnames(design_matrix)
)

# create DGEList object
experiment_dgelist <- DGEList(
  salmon_quantification$counts,
  group = metadata$condition
)

# Filter out lowly expressed genes
dim(experiment_dgelist)
keep <- filterByExpr(experiment_dgelist, design_matrix)
experiment_dgelist <- experiment_dgelist[keep,]
dim(experiment_dgelist)

# Normalise the libraries
boxplot(cpm(experiment_dgelist, log = TRUE))
experiment_dgelist <- normLibSizes(experiment_dgelist, method = "TMM")
boxplot(cpm(experiment_dgelist, log = TRUE))

plotMDS(cpm(experiment_dgelist, log = TRUE))

# Calculate mean lcpm per sample for each gene
lcpm <- cpm(experiment_dgelist, log = TRUE)

lcpm <- lcpm %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(!gene, names_to = "condition", values_to = "lcpm") %>%
  mutate(condition = experiment_dgelist$samples[condition, "group"]) %>%
  group_by(gene, condition) %>%
  summarise(lcpm = mean(lcpm))


# Remove heteroscedascity
experiment_voom <- voom(experiment_dgelist, design_matrix, plot = TRUE)
experiment_voom <- lmFit(experiment_voom, design_matrix)
experiment_voom <- contrasts.fit(experiment_voom, contrasts = contrast_matrix)
experiment_voom <- eBayes(experiment_voom)
plotSA(experiment_voom)

summary(decideTests(experiment_voom))

tfit <- treat(experiment_voom, lfc = 1)
dt <- decideTests(experiment_voom, p.value = 0.01, adjust.method = "BH")
summary(dt)

# a little bit weird, but merge all the results into one dataframe
treats <- map(
  1:ncol(dt), ~topTreat(tfit, coef = .x, n = Inf) %>%
  mutate(coef = colnames(dt)[.x]) %>%
  rownames_to_column(var = "gene")
  ) %>%
  list_rbind() %>%
  mutate(status = case_when(
    adj.P.Val < 0.01 & logFC > 0 ~ "up-regulated",
    adj.P.Val < 0.01 & logFC < 0 ~ "down-regulated",
    TRUE ~ "not significant"
  ))

treats_summary <- treats %>%
  group_by(coef, status) %>%
  summarise(n = n())

de_genes <- treats %>%
  filter(status %in% c("up-regulated", "down-regulated")) %>%
  pull(gene)

lcpm_de_genes <- lcpm %>%
  filter(gene %in% de_genes) %>%
  pivot_wider(names_from = condition, values_from = lcpm) %>%
  column_to_rownames(var = "gene")

png(file = "rnaseq_heatmap.png", width = 4, height = 10, units = "in", res = 300)
plot <- Heatmap(lcpm_de_genes)
draw(plot)
dev.off()

