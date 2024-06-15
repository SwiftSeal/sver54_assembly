import polars as pl
import logging
from rich.progress import track
from rich.logging import RichHandler

FORMAT="%(message)s"
logging.basicConfig(
    level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
)
log = logging.getLogger("rich")

def calculate_exon_mean_methylation(gene_id, gene_gff, methylation):
    # Variables for mean calculation
    sum_methylation = 0
    sum_count = 0

    exon_methylation_values = (
        gene_gff
        .filter((pl.col("gene_id") == gene_id) & (pl.col("type") == "exon"))
        .join(methylation, on="chrom")
        .filter(
            (pl.col("pos") >= pl.col("start")) & (pl.col("pos") <= pl.col("end"))
        )
        .select("methylation")
    )

    if len(exon_methylation_values) == 0:
        return None
    else:
        return exon_methylation_values.mean()

log.info("Reading gene GFF")
gene_gff = (
    pl.read_csv(
        "results/final_annotation/final_annotation.gff",
        separator="\t",
        has_header=False,
        comment_prefix="#",
    )
    .with_columns(
        pl.when(pl.col("column_3").str.contains("gene"))
        .then(pl.col("column_9").str.extract(r"ID=(.*?)$"))
        .otherwise(pl.col("column_9").str.extract(r"gene_id=(.*?);"))
        .alias("gene_id"),
    )
    .select(
        chrom=pl.col("column_1"),
        start=pl.col("column_4"),
        end=pl.col("column_5"),
        type=pl.col("column_3"),
        gene_id=pl.col("gene_id"),
    )
)

gene_methylation = {"gene": [], "type": [], "exon_mean_methylation": []}

log.info("Extracting gene IDs")
genes = (
    gene_gff
    .filter(pl.col("type") == "gene")
    .select("gene_id")
    .to_series()
    .to_list()
)

for methylation_type in ["CG", "CHG", "CHH"]:
    log.info(f"Reading {methylation_type} methylation")
    methylation = (
        pl.read_csv(
            f"results/deepsignal/freq.{methylation_type}.tsv",
            separator="\t",
            has_header=False,
        )
        .select(
            chrom=pl.col("column_1"),
            pos=pl.col("column_2"),
            methylation=pl.col("column_10"),
        )
    )

    log.info(f"Calculating mean methylation for {methylation_type}")
    for gene_id in track(genes):
        gene_methylation["gene"].append(gene_id)
        gene_methylation["type"].append(methylation_type)
        gene_methylation["exon_mean_methylation"].append(calculate_exon_mean_methylation(gene_id, gene_gff, methylation))

log.info("Saving DataFrame")
gene_methylation = pl.DataFrame(gene_methylation)
gene_methylation.write_csv("results/deepsignal/gene_exon_methylation.tsv")