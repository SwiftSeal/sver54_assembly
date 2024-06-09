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
library(topGO)

# Import data ------------------------------------------------------------------

# Experiment metadata (samples, repeats, etc)
metadata <- read.table("../config/rna_seq_samples.tsv", header = TRUE)

# Get transcript to gene relationship table from gff file
# Helixer inserted genes are mRNA, BRAKER3 or transcript
transcript_to_gene <- read.table(
  "../results/final_annotation/final_annotation.gff",
  skip = 2,
  sep = "\t"
  ) %>%
  filter(V3 == "transcript" | V3 == "mRNA") %>%
  mutate(
    gene = str_extract(V9, "Parent=([^;]+)", group = 1),
    transcript = str_extract(V9, "ID=([^;]+);", group = 1)
  ) %>%
  dplyr::select(c(transcript, gene)) %>%
  filter(!str_detect(transcript, "agat"))
  
# Salmon quantification files
sample_filepaths <- list.files("../results/salmon/")

salmon_quantification <- tximport(
  file.path("../results/salmon/", sample_filepaths, "quant.sf"),
  type = "salmon",
  tx2gene = transcript_to_gene,
  countsFromAbundance = "lengthScaledTPM"
)

# Save tpm values from salmon prior to filtering

tpm <- as.data.frame(salmon_quantification$counts)
colnames(tpm) <- sample_filepaths
tpm <- tpm %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(!gene,, names_to = "sample", values_to = "tpm") %>%
  mutate(condition = str_extract(sample, "^(.*)_Rep.*$", group = 1)) %>%
  group_by(gene, condition) %>%
  summarise(tpm = mean(tpm))
write.table(tpm, file = "../results/tpm.tsv")

# To be honest I have no idea wtf this is
# Retain GO ontology hit for highest scoring isoform
go_terms <- read.table(
  "../results/eggnog/eggnog.emapper.annotations",
  skip = 5,
  sep = "\t",
  fill = TRUE
  ) %>%
  mutate(gene = gsub(".t\\d", "", V1)) %>%
  group_by(gene) %>%
  filter(V4 == max(V4)) %>%
  dplyr::select(gene, V10)

go_terms <- as.data.frame(go_terms)

map <- go_terms[,2]
names(map) <- gsub(" ", "", go_terms[,1])
geneID2GO <- lapply(map, function(x) gsub(" ", "", strsplit(x, split = ",")[[1]]))
  

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
colnames(experiment_dgelist) <- sample_filepaths

# Filter out lowly expressed genes

dim(experiment_dgelist)
keep <- filterByExpr(experiment_dgelist, design_matrix)
experiment_dgelist <- experiment_dgelist[keep,]
dim(experiment_dgelist)

# Normalise the libraries
png(
  filename = "../results/rnaseq_counts1.png",
  width = 300,
  height = 100,
  units = "mm",
  res = 300
)
boxplot(cpm(experiment_dgelist, log = TRUE))
dev.off()
experiment_dgelist <- normLibSizes(experiment_dgelist, method = "TMM")
png(
  filename = "../results/rnaseq_counts2.png",
  width = 300,
  height = 100,
  units = "mm",
  res = 300
)
boxplot(cpm(experiment_dgelist, log = TRUE))
dev.off()

png(
  filename = "../../pandoc-thesis/figures/rnaseq_mds.png",
  width = 150,
  height = 150,
  units = "mm",
  res = 600
)
plotMDS(cpm(experiment_dgelist, log = TRUE))
dev.off()

# Calculate mean lcpm per sample for each gene
lcpm <- cpm(experiment_dgelist, log = TRUE)

lcpm <- lcpm %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(!gene, names_to = "condition", values_to = "lcpm") %>%
  mutate(condition = experiment_dgelist$samples[condition, "group"]) %>%
  group_by(gene, condition) %>%
  summarise(lcpm = mean(lcpm))

write.table(
  lcpm,
  "../results/lcpm.tsv",
  row.names = FALSE,
  sep = "\t"
)


# Remove heteroscedascity
png(
  filename = "../../pandoc-thesis/figures/rnaseq_hetero.png",
  width = 150,
  height = 100,
  units = "mm",
  res = 600
)
par(mfrow=c(1,2))
experiment_voom <- voom(experiment_dgelist, design_matrix, plot = TRUE)
experiment_voom <- lmFit(experiment_voom, design_matrix)
experiment_voom <- contrasts.fit(experiment_voom, contrasts = contrast_matrix)
experiment_voom <- eBayes(experiment_voom)
plotSA(experiment_voom)
dev.off()

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

write.table(
  treats,
  "../results/differential_expression.tsv",
  row.names = FALSE,
  sep = "\t",
)

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

png(file = "rnaseq_heatmap.png", width = 4, height = 16, units = "in", res = 600)
plot <- Heatmap(lcpm_de_genes, col = magma(100))
draw(plot)
dev.off()

# GO analysis ------------------------------------------------------------------

#geneNames <- names(geneID2GO)
#
#myInterestingGenes <- treats %>%
#  filter(coef == "Temperature25vs4") %>%
#  pull(gene)
#
#
#geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
#names(geneList) <- geneNames
#str(geneList)
#
#GOdata <- new(
#  "topGOdata",
#  ontology = "MF",
#  allGenes = geneList,
#  gene2GO = geneID2GO,
#  annot = annFUN.gene2GO,
#  nodeSize = 5
#  )
#
#resultFisher <- runTest(GOdata, algorithim = "classic", statistic = "fisher")
#
#allRes <- GenTable(
#  GOdata,
#  classicFisher = resultFisher
#)
#View(allRes)
