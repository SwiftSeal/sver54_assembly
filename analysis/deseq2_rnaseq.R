library(DESeq2)
library(tximport)
library(readr)
library(stringr)
library(dplyr)

# Import data ------------------------------------------------------------------

tx2gene <- read_tsv("results/nf_rnaseq/tx2gene.tsv", col_names = FALSE)

# Salmon quantification files
sample_filepaths <- fs::path("results", "salmon") |>
  list.files()

samples <- read_csv("config/rnaseq.csv") |>
  mutate(condition = str_remove(sample, "_REP\\d"))

sample_filepaths <- fs::path("results/nf_rnaseq/", samples$sample, "quant.sf")

txi <- tximport(
  sample_filepaths,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)


ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)


vsd <- vst(dds)

plotPCA(vsd)

dds <- DESeq(ddsTxi, betaPrior = TRUE)
res <- results(
  dds,
  contrast = c("condition", "TEMPERATURE_4", "TEMPERATURE_25"),
)
summary()

plotMA(res.shrink)
