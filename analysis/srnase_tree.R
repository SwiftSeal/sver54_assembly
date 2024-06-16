library(ggtree)
library(gggenomes)
library(dplyr)
library(stringr)
library(tidyr)

# Plot tree of refseq hits -----------------------------------------------------

# Import the tree - then zoom into the clade local to candidate
tree <- read.tree("../results/srnase_refseq_hits.mafft.fa.raxml.bestTree")
zoomed <- tidytree::tree_subset(tree, "g3195.t1", levels_back = 9)

# Extract useful data from the fasta headers
headers <- readLines("../results/srnase_refseq_hits.fa")
headers <- headers[grep(">", headers)]
ids <- str_extract(headers, ">(.*?) (.*?) \\[(.*?)\\]", group = 1)
descriptions <- str_extract(headers, ">(.*?) (.*?) \\[(.*?)\\]", group = 2)
species <- str_extract(headers, ">(.*?) (.*?) \\[(.*?)\\]", group = 3)

# Make a table of this data, and create a pretty label of tips
metadata <- tibble(ids, descriptions, species) %>%
  mutate(
    formatted_label = str_c(
      ids,
      "~italic(",
      str_replace(species, " ", "~"),
      ")"
      )
    )

# Replace candidate name with something more explicit
metadata <- metadata %>%
  mutate(
    ids = replace_na(ids, "g3195.t1"),
    formatted_label = replace_na(formatted_label, "My~gene")
  )

tree_plot <- ggtree(zoomed) %<+% metadata +
  geom_tiplab(align = TRUE, parse = TRUE, size = 3, aes(label = formatted_label)) +
  xlim(0, 4) +
  geom_treescale(x = 0, y = 20, width = 0.8)


ggsave("srnase_tree.png", tree_plot, width = 5.9, height = 4, dpi = 600)

# Plot local gene environment --------------------------------------------------

genes <- read_gff3(
  "../results/final_annotation/final_annotation.gff"
)

earlgrey <- read.table(
  "../results/earlgrey/solanum_verrucosum_summaryFiles/solanum_verrucosum.filteredRepeats.gff"
) %>%
  select(seq_id=V1, start=V4, end=V5, type = V3)

edta <- read.table(
  "../results/edta/final_assembly.fa.mod.EDTA.TEanno.gff3"
) %>%
  select(seq_id=V1, start=V4, end=V5, type = V3)

locus <- genes %>%
  filter(feat_id == "g3195.t1") %>%
  mutate(
    start = start - 5000,
    end = end + 5000,
    length = end + 10000
  ) %>%
  select(seq_id, length, start, end) %>%
  distinct()

gggenomes(seqs = locus, genes = genes %>% filter(type != "mRNA") %>% mutate(type = ifelse(type == "transcript", "mRNA", type)), feats = edta) +
  geom_seq() +
  geom_gene(intron_types = "mRNA") +
  geom_feat(aes(colour = type))
