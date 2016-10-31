"""
Do Tiling for the 2-D MESA UQ project.
"""
import numpy as np
from Tiling import Point, Tile, Domain

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
dom = Domain(pointlist, lo, hi)

# Tile Domain
dom.do_domain_tiling()

# List Tiles and Bounds
for i, t in enumerate(dom.tiles):
    print('Tile {}:'.format(i))
    print('-- lo: {}'.format(t.lo))
    print('-- hi: {}'.format(t.hi))
