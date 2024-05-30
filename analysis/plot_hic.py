import cooler
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan

clr = cooler.Cooler('results/hic/contact_maps/cool/sample1.mcool::/resolutions/320000')

chromstarts = []
for i in clr.chromnames:
    chromstarts.append(clr.extent(i)[0])

cg = adaptive_coarsegrain(clr.matrix(balance=True)[:], clr.matrix(balance=False)[:], cutoff = 3)


cgi = interp_nan(cg)

f, axs = plt.subplots(
    figsize=(10, 10)
)

norm=LogNorm(vmax=0.05)
im = axs.matshow(cgi, norm = norm, cmap = 'magma')
axs.set(xticks = chromstarts, xticklabels = clr.chromnames)



# save
plt.savefig('sample1.png')