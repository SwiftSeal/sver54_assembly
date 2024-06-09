make_tracks_file \
    --trackFiles results/cenh3/aln.bw \
    results/deepsignal/freq.CG.bw \
    results/deepsignal/freq.CHG.bw \
    results/deepsignal/freq.CHH.bw \
    results/TRASH/TRASH_final_assembly.sorted.small.bed \
    results/earlgrey/solanum_verrucosum_EarlGrey/solanum_verrucosum_summaryFiles/LTR.sorted.bed \
    -o results/genometrack.ini

pyGenomeTracks --tracks results/genometrack.ini --BED results/cenh3/centromeres.bed --outFileName results/centromeres.png