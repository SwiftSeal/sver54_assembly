multiBigwigSummary BED-file \
  -b \
    "/mnt/shared/scratch/msmith/solanum_verrucosum/results/deepsignal/freq.CG.bw" \
    "/mnt/shared/scratch/msmith/solanum_verrucosum/results/deepsignal/freq.CHG.bw" \
    "/mnt/shared/scratch/msmith/solanum_verrucosum/results/deepsignal/freq.CHH.bw" \
  --BED "~/scratch/solanum_verrucosum/results/genes/solanum_verrucosum.gtf" \
  --metagene
  -out results/gene_methylation.npz \
  --outRaw results/gene_methylation.txt
