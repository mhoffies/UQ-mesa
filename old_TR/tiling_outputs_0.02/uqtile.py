"""
Do Tiling for the 2-D MESA UQ project.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle
from Tiling import Point, Tile, Domain, DMCycle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str, help='Name of the input csv file containing (x, y, value) data series.')
args = parser.parse_args()

# Read Data
raw_data = np.genfromtxt(args.infile, delimiter=',', skip_header=1)
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
dom.do_domain_tiling(gnr_thresh=0.02, attempt_virtual_shrink=True)
