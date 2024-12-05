import cooler
import cooltools
import matplotlib.pyplot as plt
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan
from matplotlib.colors import LogNorm

# Load the Hi-C data
#cooler.fileops.list_coolers('results/hic/contact_maps/cool/sample1.mcool')
clr = cooler.Cooler(
    'results/hic/contact_maps/cool/sample1.mcool::resolutions/640000'
)

cg = adaptive_coarsegrain(
    clr.matrix(balance=True)[:],
    clr.matrix(balance=False)[:],
    cutoff=3,
    max_levels=8
)

cgi = interp_nan(cg)

norm = LogNorm(vmax=0.01)

chromstarts = []
for i in clr.chromnames:
    print(f'{i} : {clr.extent(i)}')
    chromstarts.append(clr.extent(i)[0])

f, ax = plt.subplots(
    figsize=(8, 8),
    dpi=600,
)

im = ax.matshow(cgi, norm=norm)

ax.set(xticks=chromstarts, xticklabels=clr.chromnames,
       xlabel='position, chrom#', ylabel='position, bin#')
ax.xaxis.set_label_position('top')
f.colorbar(im, fraction=0.046, pad=0.04)
plt.savefig('results/hic_plot.png')