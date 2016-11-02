"""
Do Tiling for the 2-D MESA UQ project.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle
from Tiling import Point, Tile, Domain, DMCycle

# Read Data
dfname = 'output2.csv'
raw_data = np.genfromtxt(dfname, delimiter=',', skip_header=1)
# Each element of data is a row from the csv file, so convert to columns
data = np.transpose(raw_data)
# data[0] = Blocker factors
xvec = data[0]
# data[1] = Reimers factors
yvec = data[1]
# data[2] = CO WD Mass
zvec = data[2]

# Create list of Points
pointlist = []
for x, y, z in zip(xvec, yvec, zvec):
    p = Point([x,y], z)
    pointlist.append(p)

# Get bounds on the x,y domain
lo = [np.amin(xvec), np.amin(yvec)]
hi = [np.amax(xvec), np.amax(yvec)]

# Form Domain
dom = Domain(points=pointlist, lo=lo, hi=hi)

# Tile Domain
dom.do_domain_tiling(gnr_thresh=0.05)

# Plot Domain
fig = plt.figure()
ax = fig.add_subplot(111)
# Plot Tile outlines
linestyle_options = ['-', '--', '-.', ':']
ls_cycler = DMCycle(len(linestyle_options))
for t in dom.tiles:
    ax.add_patch(Rectangle((t.lo[0], t.lo[1]), t.hi[0]-t.lo[0], t.hi[1]-t.lo[1], facecolor='None',
                           edgecolor='orange',
                           linewidth=1.5,
                           linestyle=linestyle_options[ls_cycler.cycle()]))
# Plot points inside Tiles
points_x = []
points_y = []
points_v = []
for j, t in enumerate(dom.tiles):
    for i, p in enumerate(t.points):
        points_x.append(p.r[0])
        points_y.append(p.r[1])
        points_v.append(p.v)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
cmap = mpl.cm.viridis
bounds = [0.500,0.525,0.550,0.575,0.600]
norm = mpl.colors.Normalize(vmin=np.amin(bounds),vmax=np.amax(bounds))
img = ax.scatter(points_x, points_y, c=points_v, cmap=cmap, norm=norm)
plt.colorbar(img, cmap=cmap, cax=cax, ticks=bounds, norm=norm, label='$M_{WD}~(M_{\odot})$')
# Plot points outside Tiles
for i, p in enumerate(dom.scratch_points):
    ax.scatter(p.r[0], p.r[1], color='red')
ax.set_ylabel('$\eta_R$')
ax.set_xlabel('$\eta_B$')
plt.tight_layout()
plt.savefig('tiled_domain.eps')
