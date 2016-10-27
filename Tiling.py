"""
Tile a 2-D domain containing Point objects depending on a decision function.

Donald Willcox
"""
import numpy as np
from Plane_nd import Plane_nd

class Point(object):
    def __init__(self, r=[], v=None):
        # Position in n-D space
        self.r = np.array(r)
        # Value of point (Scalar)
        self.v = v
        # n-D Mask: Boundary Edge represented by Point
        self.bedg = [None for i in range(len(r))]
        # n-D Mask: Boundary Type represented by Point
        self.bdir = self.bcon

    def dist_to_pt(self, b):
        # Get distance between this Point and Point b
        dr = np.array(self.r) - np.array(b.r)
        return np.sqrt(np.sum(dr**2))


class Plane(object):
    def __init__(self, points=None, dm=None):
        self.cpars = None # Length n+1 for n-D space
        self.dm = dm
        self.resd  = None
        self.norm_resd = None
        self.geom_norm_resd = None
        self.compute_pars(points)
        
    def compute_pars(self, points):
        if not points:
            return
        ivars = np.array([p.r for p in points])
        dvars = np.array([p.v for p in points])
        fitter = Plane_nd(ivars, dvars, self.dm)
        popt, pcov = fitter.dolsq()
        dpfit = fitter.fplane(ivars, popt)
        self.resd = dvars - dpfit
        self.norm_resd = self.resd/dvars
        self.geom_norm_resd = np.sqrt(np.sum(self.norm_resd**2))
        
class Tile(objects):
    def __init__(self, points=[], lo=[], hi=[], dm=None):
        self.children = []
        self.points = points
        self.lo = lo
        self.hi = hi
        self.dm = dm
        self.nres_threshold = 0.5
        self.min_contain_points = 3

    # # 2D
    # def get_vertices(self):
    #     # Return the vertices of the rectangle defined by lo, hi
    #     vertices = []
    #     vertices.append(self.lo)
    #     vertices.append([self.lo[0], self.hi[1]])
    #     vertices.append(self.hi)
    #     vertices.append([self.hi[0], self.lo[1]])
    #     return vertices
    
    def get_enclosed_points(self, lo, hi):
        # Return list of self.points within [lo, hi]
        inpts = []
        for pt in self.points:
            pt_in = True
            for di in range(dm):
                if pt.r[di] < lo[di] or pt.r[di] > hi[di]:
                    pt_in = False
                    break
            if pt_in:
                inpts.append(pt)
        return inpts

    def get_subtile(self, lo, hi):
        # Return a Tile object corresponding to a subtile of self.
        # Return None if Error.
        
        # First, check domain partitioning
        if len(lo) != self.dm:
            return None
        if len(hi) != self.dm:
            return None
        for di in range(self.dm):
            if lo[di] < self.lo[di] or lo[di] > self.hi[di]:
                return None
            if hi[di] < self.lo[di] or hi[di] > self.hi[di]:
                return None
        # Now get the points within [lo, hi]
        inpts = self.get_enclosed_points(lo, hi)
        # Create and return sub-Tile
        stile = Tile(inpts, lo, hi, self.dm)
        return stile

    def get_geom_norm_resd(self):
        # Returns geometric mean of normalized residuals
        # between the points in the Tile and a Plane fit.
        p = Plane(self.points, self.dm)
        return p.geom_norm_resd

    # # 2D
    # def determine_subdivide(self):
    #     if len(self.points) == 3:
    #         return False
    #     elif len(self.points) < 3:
    #         print('ERROR: Tile contains less than 3 points')
    #         exit()
    #     geom_norm_resd = self.get_geom_norm_resd()
    #     if geom_norm_resd > self.nres_threshold:
    #         return True
    #     else:
    #         return False

class Domain(object):
    def __init__(self, tiles=[], points=[], lo=[], hi=[]):
        # The Domain is just a set of Point objects
        # and functions for tiling them into a set of Tile objects.
        
        self.tiles = tiles
        self.lo = lo
        self.hi = hi
        self.dm = None

        if lo and hi and len(lo) != len(hi):
            print('ERROR: lo and hi supplied with incongruous dimensions.')
            exit()
        else:
            self.dm = len(lo)
        
        # Create a tile if only points, lo and hi are given
        if lo and hi and points and not tiles:
            t = Tile(points, lo, hi, self.dm)
            self.tiles.append(t)

    def get_distal_point(self, refpt, points):
        # Get the most distal point from refpt among points
        dmax = 0.0
        pdst = None
        for p in points:
            dp = self.get_dpt(refpt, p)
            if dp > dmax:
                dmax = dp
                pdst = p
        return pdst

    def get_nearest_point(self):
        # Get the point nearest to a side of an n-D rectangle

    def extend_rectangle(self, vert):
        for cpt in closure:
            uplim = self.get_up_limit_rectangle(cpt)
            inpts = self.get_enclosed_points(self.points, cpt, uplim)

    def partition(self):
        redo = True
        while(redo):
            for t in self.tiles:
                while t.determine_subdivide():
                    lo = t.lo
                    hi = t.hi
                    divs = []
                    for di in range(t.dm):
                        dvi = 0.5*(lo[di] + hi[di])
                        lo_dn = lo
                        hi_dn = hi
                        hi_dn[di] = dvi
                        lo_up = lo
                        hi_up = hi
                        lo_up[di] = dvi
                        tile_dn = t.get_subtile(lo_dn, hi_dn)
                        tile_up = t.get_subtile(lo_up, hi_up)
                        gnrd = tile_dn.get_geom_norm_resd() + tile_up.get_geom_norm_resd()
                        dd = {'gnrd': gnrd, 'tdn': tile_dn, 'tup': tile_up}
                        divs.append(dd)
                    # Figure out which subdivision did best.
                    

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
    p = Point(x,y,z)
    pointlist.append(p)

# Get bounds on the x,y domain
lo = [np.amin(xvec), np.amin(yvec)]
hi = [np.amax(xvec), np.amax(yvec)]

# Form Domain
dom = Domain(pointlist, lo, hi)

        
