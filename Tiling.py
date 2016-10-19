"""
Tile a 2-D domain containing Point objects depending on a decision function.

Donald Willcox
"""
import numpy as np
from scipy.optimize import curve_fit

class Point(object):
    def __init__(self, r=[], v=None)
        self.r = r # Position in n-D space
        self.v = v # Value of point

class Plane(object):
    def __init__(self, points=None):
        self.c0 = None
        self.cx = None
        self.cy = None
        self.resd = None
        self.norm_resd = None
        self.geom_norm_resd = None
        self.compute_pars(points)
        
    def fplane(self, xy, cy, cx, c0):
        x = xy[0]
        y = xy[1]
        return c0 + cx * x + cy * y
    
    def compute_pars(self, points):
        if not points:
            return
        xvec = []
        yvec = []
        zvec = []
        for p in points:
            xvec.append(p.x)
            yvec.append(p.y)
            zvec.append(p.z)
        xvec = np.array(xvec)
        yvec = np.array(yvec)
        zvec = np.array(zvec)
        ivars = np.array([xvec, yvec])
        popt, pcov = curve_fit(self.fplane, ivars, zvec)
        pstd = np.sqrt(np.diag(pcov))
        self.cy, self.cx, self.c0 = popt
        zfit = self.fplane(ivars, self.cy, self.cx, self.c0)
        self.resd = zvec - zfit
        self.norm_resd = self.resd/zvec
        self.geom_norm_resd = np.sqrt(np.sum(self.norm_resd**2))
        
class Tile(objects):
    def __init__(self, points=[], lo=[0,0], hi=[1,1]):
        self.children = []
        self.points = points
        self.lo = lo
        self.hi = hi
        self.dm = 2 # Number of dimensions
        self.nres_threshold = 0.5
        self.min_contain_points = 3
        
    def get_vertices(self):
        # Return the vertices of the rectangle defined by lo, hi
        vertices = []
        vertices.append(self.lo)
        vertices.append([self.lo[0], self.hi[1]])
        vertices.append(self.hi)
        vertices.append([self.hi[0], self.lo[1]])
        return vertices
    
    def get_enclosed_points(self, lo, hi):
        # Return list of self.points within [lo, hi]
        inpts = []
        for pt in self.points:
            pt_in = True
            for di in xrange(dm):
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
        for di in xrange(self.dm):
            if lo[di] < self.lo[di] or lo[di] > self.hi[di]:
                return None
            if hi[di] < self.lo[di] or hi[di] > self.hi[di]:
                return None
        # Now get the points within [lo, hi]
        inpts = self.get_enclosed_points(lo, hi)
        # Create and return sub-Tile
        stile = Tile(inpts, lo, hi)
        return stile

    def get_dpt(self, a, b):
        dr = np.array(a.r) - np.array(b.r)
        return np.sqrt(np.sum(dr**2))

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

    def get_distal_x(self, refpt, points):
        # Get most distal point in x direction
        dmax = 0.0
        pdst = None
        for p in points:
            dp = abs(p.x - refpt.x)
            if dp > dmax:
                dmax = dp
                pdst = p
        return pdst

    def get_distal_y(self, refpt, points):
        # Get most distal point in y direction
        dmax = 0.0
        pdst = None
        for p in points:
            dp = abs(p.y - refpt.y)
            if dp > dmax:
                dmax = dp
                pdst = p
        return pdst
            
    def extend_rectangle(self, vert):
        for cpt in closure:
            uplim = self.get_up_limit_rectangle(cpt)
            inpts = self.get_enclosed_points(self.points, cpt, uplim)
            
    def determine_subdivide(self):
        if len(self.points) == 3:
            return False
        elif len(self.points) < 3:
            print('ERROR: Tile contains less than 3 points')
            exit()
        p = Plane(self.points)
        if p.geom_norm_resd > self.nres_threshold:
            return True
        else:
            return False

class Domain(object):
    def __init__(self, tiles=[], points=[], lo=[], hi=[]):
        # The Domain is just a set of Tile objects
        # and functions for managing them.
        
        self.tiles = tiles
        self.lo = lo
        self.hi = hi
        
        # Create a tile if only points, lo and hi are given
        if lo and hi and points and not tiles:
            t = Tile(points, lo, hi)
            self.tiles.append(t)

    def partition(self):
        redo = True
        while(redo):
            for t in self.tiles:
                while t.determine_subdivide():
                    lo = t.lo
                    hi = t.hi
                    for di in xrange(t.dm):
                        dvi = 0.5*(lo[di] + hi[di])
                        lo_dn = lo
                        hi_dn = hi
                        hi_dn[di] = dvi
                        lo_up = lo
                        hi_up = hi
                        lo_up[di] = dvi
                        tile_dn = t.get_subtile(lo_dn, hi_dn)
                        tile_up = t.get_subtile(lo_up, hi_up)

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

        
