library(tidyverse)

# Plot the S-RNase gene with upstream TEs and methylation

results_path = "Z://scratch/sver54_assembly/results"

srnase_id <- "cpc54_3362.1"

genes <- read_tsv(
  file.path(results_path, "annotation/cpc54.gene.gff"),
  col_names = FALSE,
  comment = "#"
) |>
  filter(str_detect(X9, srnase_id))

srnase_start <- 47986461
srnase_end <- 47985683

earlgrey <- read_tsv(
  file.path(results_path, "hite/HiTE.gff"),
  col_names = FALSE,
  comment = "#"
) |>
  filter(X4 > srnase_start, X5 < srnase_start + 2000, X1 == "chr01")

cg <- read_tsv(
  file.path(results_path, "methylation/cpc54.cg.bedgraph"),
  col_names = FALSE,
) |>
  filter(X2 > srnase_end, X3 < srnase_start + 2000, X1 == "chr01")

chg <- read_tsv(
  file.path(results_path, "methylation/cpc54.chg.bedgraph"),
  col_names = FALSE,
) |>
  filter(X2 > srnase_end, X3 < srnase_start + 2000, X1 == "chr01")

chh <- read_tsv(
  file.path(results_path, "methylation/cpc54.chh.bedgraph"),
  col_names = FALSE,
) |>
  filter(X2 > srnase_end, X3 < srnase_start + 2000, X1 == "chr01")

srnase_plot <- ggplot() +
  geom_rect(data = genes, aes(xmin = X4, xmax = X5, ymin = 0, ymax = -.1, fill = X3)) +
  geom_rect(data = earlgrey, aes(xmin = X4, xmax = X5, ymin = 0, ymax = -.1, fill = X9)) +
  geom_rect(data = cg, aes(xmin = X2, xmax = X3, ymin = 0, ymax = 0 + X4), fill = "#4E79A7") +
  geom_rect(data = chg, aes(xmin = X2, xmax = X3, ymin = 0, ymax = 0 + X4), fill = "#F28E2B") +
  geom_rect(data = chh, aes(xmin = X2, xmax = X3, ymin = 0, ymax = 0 + X4), fill = "#E15759") +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank()) +
  labs(x = "Genome position (bp)", y = "Methylation") +
  theme(legend.position = "none")


ggsave("plots/srnase_plot.pdf", srnase_plot, width = 5.9, height = 1)
```

```{r}
srnase_aa <- readAAStringSet("../results/OG0023638.fa")
msaPrettyPrint(msa(srnase_aa))

slf_aa <- readAAStringSet("../results/OG0000138.fa")
msa(slf_aa)
```

```{r message=FALSE}
get_upstream <- function(bed, edta, gene, boundary) {
  helixer <- read_delim(bed, col_names = FALSE)
  
  start <- as.numeric(helixer[helixer$X4 == gene, "X2"])
  
  df <- read_tsv(edta, col_names = FALSE, comment = "#") |>
    filter(X1 == "chr01", X4 > start - boundary, X4 < start) |>
    dplyr::select(xmin = X4, xmax = X5, feature = X3) |>
    mutate(
      gene = gene,
      xmin = xmin - start,
      xmax = xmax - start
    )
  
  return(df)
}

dfs <- list()

for (gene in srnase$tip.label) {
  genome <- str_extract(gene, "(.*)_chr.*$", group = 1)
  
  df <- get_upstream(
    paste0("../results/pangenomics_annot/", genome, ".bed"),
    paste0("../results/pangenomics_annot/", genome, ".fa.mod.EDTA.TEanno.gff3"),
    gene,
    5000
  )
  
  dfs[[length(dfs) + 1]] <- df
}

df <- bind_rows(dfs)

df <- df |>
  mutate(
    feature = case_when(
      str_detect(feature, "LTR") ~ "LTR",
      str_detect(feature, "TIR") ~ "TIR",
      .default = NA
    )
  ) |>
  drop_na()
```


```{r}
earlgrey <- read_tsv(
  "../results/earlgrey/solanum_verrucosum_summaryFiles/solanum_verrucosum.filteredRepeats.gff",
  col_names = FALSE
) |>
  filter(X1 == "chr01", X4 > 47986461, X5 < 47986461 + 3000) |>
  mutate(X4 = X4 - 47986461, X5 = X5 - 47986461)

edta <- read_tsv(
  "../results/edta/final_assembly.fa.mod.EDTA.TEanno.gff3",
  col_names = FALSE,
  comment = "#"
) |>
  filter(X1 == "chr01", X4 > 47986461, X5 < 47986461 + 3000) |>
  mutate(X4 = X4 - 47986461, X5 = X5 - 47986461)

annotation_plot = ggplot() +
  geom_segment(data = earlgrey, x = 0, xend = 3000, y = "Earlgrey", yend = "Earlgrey", colour = "#333333", size = 1) +
  geom_segment(data = earlgrey, x = 0, xend = 3000, y = "EDTA", yend = "EDTA", colour = "#333333", size = 1) +
  geom_segment(data = earlgrey, aes(x = X4, xend = X5, y = "Earlgrey", yend = "Earlgrey", colour = X3), size = 5) +
  geom_segment(data = edta, aes(x = X4, xend = X5, y = "EDTA", yend = "EDTA", colour = X3), size = 5) +
  scale_colour_tableau(name = "Annotation") +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank()) +
  theme(legend.position = "top") +
  labs(y = NULL, x = "S-RNase upstream (bp)")

ggsave("plots/srnase_upstream.png", width = 5.9, height = 1.5, dpi = 600)
  
  
```

```{python}

```

```{r}
plot <-  df |>
  drop_na() |>
  mutate(
    feature = case_when(
      str_detect(feature, "LTR") ~ "LTR",
      str_detect(feature, "TIR") ~ "TIR",
      .default = "Other"
    )
  ) |>
  ggplot(aes(xmin = xmin, xmax = xmax, fill = feature, y = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ gene, ncol = 1, scales = "free")

ggsave("plots/srnase.png", plot, width = 10, height = 10)
```

```{r}
gene_gff <- read_tsv(
  "../results/annotation/cpc54.gene.gff",
  col_names = FALSE,
  comment = "#"
) |>
  filter(str_detect(X9, "cpc54_3362.1"))

start = 47985683
end = 47986461

te_gff <- read_tsv(
  "../results/annotation/"
)

cg <- read_tsv("../results/methylation/freq.CG.bedgraph", col_names = NA) |>
  filter(X1 == "chr01", X2 > start - 5000, X3 < end + 5000)

chg <- read_tsv("../results/methylation/freq.CHG.bedgraph", col_names = NA) |>
  filter(X1 == "chr01", X2 > start - 5000, X3 < end + 5000)

chh <- read_tsv("../results/methylation/freq.CHH.bedgraph", col_names = NA) |>
  filter(X1 == "chr01", X2 > start - 5000, X3 < end + 5000)
```

```{r}
ggplot() +
  geom_rect(data = gene_gff, aes(xmin = X4, xmax = X5, ymin = 0, ymax = 1, fill = X3)) +
  geom_segment(data = cg, aes(x = X2, xend = X2, y = 1, yend = 1 + X4)) +
  geom_segment(data = chg, aes(x = X2, xend = X2, y = 2, yend = 2 + X4)) +
  geom_segment(data = chh, aes(x = X2, xend = X2, y = 3, yend = 3 + X4)) 
```

