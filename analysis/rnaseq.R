library(dplyr)
library(stringr)
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
  levels = colnames(design)
)

# create DGEList object
experiment_dgelist <- DGEList(
  salmon_quantification$counts,
  group = metadata$condition
)

# Filter out lowly expressed genes
dim(experiment_dgelist)
keep <- filterByExpr(experiment_dgelist, design)
experiment_dgelist <- experiment_dgelist[keep,]
dim(experiment_dgelist)

# Normalise the libraries
boxplot(cpm(experiment_dgelist, log = TRUE))
experiment_dgelist <- normLibSizes(experiment_dgelist, method = "TMM")
boxplot(cpm(experiment_dgelist, log = TRUE))

plotMDS(cpm(experiment_dgelist, log = TRUE))

# Remove heteroscedascity
experiment_voom <- voom(experiment_dgelist, design_matrix, plot = TRUE)
experiment_voom <- lmFit(experiment_voom, design_matrix)
experiment_voom <- contrasts.fit(experiment_voom, contrasts = contrast_matrix)
experiment_voom <- eBayes(experiment_voom)
plotSA(experiment_voom)

summary(decideTests(experiment_voom))

tfit <- treat(vfit, lfc = 1)
dt <- decideTests(tfit, p.value = 0.01, adjust.method = "BH")
summary(dt)

# a little bit weird, but merge all the results into one dataframe
treats <- map(
  1:ncol(dt),
  ~topTreat(tfit, coef = .x, n = Inf) %>%
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

# plot number of significant differential gene expressions
plot_significant <- ggplot(treats_summary %>% filter(status != "not significant"), aes(y = coef, x = n, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("up-regulated" = colours[1], "down-regulated" = colours[2]))

ggsave(plot_significant, filename = "results/significant_genes.png", width = 6, height = 3, units = "in", dpi = 300)

# plot mean expression vs logFC
plot_MA <- ggplot(treats, aes(x = logFC, y = AveExpr, colour = status)) +
  geom_point(data = treats %>% filter(status == "not significant")) +
  geom_point(data = treats %>% filter(status == "up-regulated")) +
  geom_point(data = treats %>% filter(status == "down-regulated")) +
  facet_wrap(~coef, ncol = 2) +
  theme(legend.position = c(0.75, 0.12), legend.justification = c(0.5, 0.5)) +
  scale_colour_manual(values = c("not significant" = colours[7], "up-regulated" = colours[1], "down-regulated" = colours[2]))

ggsave(plot_MA, filename = "results/MA_plot.png", width = 6, height = 6, units = "in", dpi = 300)

write.table(treats, file = "results/differential_expression.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

