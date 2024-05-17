library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)
library(stringr)
library(Biostrings)
library(ggtree)

nice_colours = c("#2596be", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

# Load and parse EarlGrey data -------------------------------------------------

tesorter_earlgrey <- read.table(
  "../results/tesorter/earlgrey.cls.tsv",
  header = TRUE,
  comment.char = "",
  sep = "\t") %>%
  mutate(X.TE = tolower(X.TE))

fasta_earlgrey <- readDNAStringSet(
  "../results/earlgrey/solanum_verrucosum_summaryFiles/solanum_verrucosum-families.fa.strained")

statistics_earlgrey<- data.frame(
  Sequence = trimws(names(fasta_earlgrey)),
  Length = width(fasta_earlgrey)) %>%
  mutate(Name = str_extract(Sequence, "^[^#/]+"),
         Classification = str_extract(Sequence, "(?<=#)[^/]+"),
         Subclass = str_extract(Sequence, "(?<=/).+")) %>%
  mutate(Classification = case_when(
    Subclass == "Helitron" ~ "Helitron",
    Classification == "rRNA" ~ "Other",
    Classification == "Satellite" ~ "Other",
    Classification == "Simple_repeat" ~ "Other",
    Classification == "snRNA" ~ "Other",
    Classification == "tRNA" ~ "Other",
    Classification == "Unknown" ~ "Other",
    .default = Classification)) %>%
  mutate(Source = "EarlGrey") %>%
  mutate(Sequence = tolower(Sequence))

coverage_earlgrey <- read.table(
  "../results/earlgrey/solanum_verrucosum_summaryFiles/solanum_verrucosum.familyLevelCount.txt",
  header = TRUE,
  comment.char = "",
  sep = "\t") %>%
  select(Sequence = name, Coverage = coverage) %>%
  mutate(Sequence = tolower(Sequence))

earlgrey <- statistics_earlgrey %>%
  left_join(tesorter_earlgrey, by = join_by("Sequence" == "X.TE")) %>%
  left_join(coverage_earlgrey, by = join_by("Sequence"))

# Load and parse EDTA data -----------------------------------------------------

tesorter_edta <- read.table(
  "../results/tesorter/tesorter.cls.tsv",
  header = TRUE,
  comment.char = "",
  sep = "\t")

fasta_edta <- readDNAStringSet(
  "../results/edta/final_assembly.fa.mod.EDTA.TElib.fa")

statistics_edta <- data.frame(
  Sequence = trimws(names(fasta_edta)),
  Length = width(fasta_edta)) %>%
  mutate(Name = str_extract(Sequence, "^[^#/]+"),
         Classification = str_extract(Sequence, "(?<=#)[^/]+"),
         Subclass = str_extract(Sequence, "(?<=/).+")) %>%
  mutate(Classification = case_when(
    Subclass == "Helitron" ~ "Helitron",
    Classification == "MITE" ~ "Other",
    .default = Classification)) %>%
  mutate(Source = "EDTA")

# EDTA has silly formatting - have to specify table cutoff
coverage_edta <- read.table(
  "../results/edta/final_assembly.fa.mod.EDTA.TEanno.sum",
  skip = 40,
  nrows = 5951) %>%
  select(Name = V1, Coverage = V3)

edta <- statistics_edta %>%
  left_join(tesorter_edta, by = join_by("Sequence" == "X.TE")) %>%
  left_join(coverage_edta, by = join_by("Name"))


# Merge ------------------------------------------------------------------------

merged <- bind_rows(earlgrey, edta)

# Plot summary statistics ------------------------------------------------------

coverage_barplot <- merged %>%
  group_by(Source, Classification) %>%
  summarise(Coverage = sum(Coverage, na.rm = TRUE)) %>%
  ggplot(aes(x = Classification, y = Coverage, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Coverage (bp)") +
  scale_y_continuous(label = label_number(scale_cut = cut_short_scale(), suffix = "bp")) +
  theme(legend.position = c(0,0)) +
  scale_fill_manual(values = nice_colours)

length_boxplot <- merged %>%
  group_by(Source, Classification) %>%
  ggplot(aes(x = Classification, y = Length, fill = Source)) +
  geom_boxplot() +
  labs(y = "Length (bp)") +
  scale_y_continuous(label = label_number(scale_cut = cut_short_scale(), suffix = "bp")) +
  theme(legend.position = "none") +
  scale_fill_manual(values = nice_colours)

count_barplot <- merged %>%
  group_by(Source, Classification) %>%
  ggplot(aes(x = Classification, fill = Source)) +
  geom_bar(position = "dodge") +
  labs(y = "Families (#)") +
  theme(legend.position = "none") +
  scale_fill_manual(values = nice_colours)

summary_plot <- plot_grid(
  coverage_barplot,
  length_boxplot,
  count_barplot,
  align = "hv",
  ncol = 1)

summary_plot

# Some numbers for the paper ---------------------------------------------------

# How many TEs could not be classified by TEsorter?
merged %>%
  group_by(Source, Classification) %>%
  summarise(value = 100*sum(is.na(Order)/n()))

# How many TEs were considered "Complete"?
merged %>%
  group_by(Source, Classification) %>%
  summarise(
    value = sum(Complete == "yes", na.rm = TRUE),
    prop = 100*sum(Complete == "yes", na.rm = TRUE)/n())

# Plot phylogenetic tree of LTRs -----------------------------------------------

tree <- read.tree("../results/tesorter/earlgrey.cls.pep.RT.aln.parstree")

# Fix tip labels so they fit form 'rnd-1_family-5'

tree$tip.label <- str_extract(tree$tip.label, "rnd-\\d+_family-\\d+")

earlgrey_named <- earlgrey %>%
  select(Name:Coverage)

earlgrey_named %>%
  group_by(Clade) %>%
  summarise(total = sum(Coverage, na.rm = TRUE)) %>% View()

tree_plot <- ggtree(tree, layout = "daylight") %<+% earlgrey_named
tree_plot %>%
  groupClade(.node = 21) + aes(colour = group)
tree_plot <- tree_plot +
  geom_text(aes(label=node)) +
  geom_tippoint(aes(colour = Superfamily))
ggsave("tree_plot.pdf", plot = tree_plot, width = 20, height = 20, units = "in")
