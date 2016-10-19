"""
Tile a 2-D domain containing Point objects depending on a decision function.

Donald Willcox
"""
import numpy as np
from scipy.optimize import curve_fit

class Point(object):
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

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
    def __init__(self, points=[], lo=[0,0], hi=[1,1], nres_threshold=0.5):
        self.children = []
        self.closure = [lo] # List of points defining lower boundary
        self.points = points
        self.lo = lo
        self.hi = hi
        self.nres_threshold = nres_threshold
        self.min_contain_points = 3
        
    def get_vertices(self):
        # Return the vertices of the rectangle defined by lo, hi
        vertices = []
        vertices.append(self.lo)
        vertices.append([self.lo[0], self.hi[1]])
        vertices.append(self.hi)
        vertices.append([self.hi[0], self.lo[1]])
        return vertices
    
    def get_up_limit_rectangle(self, vert):
        # Extend a rectangle starting at a point defining its lower vertex
        vertx = vert[0]
        verty = vert[1]
        hixlm = [self.hi[0]]
        hiylm = [self.hi[1]]
        for ctile in self.children:
            lox = ctile.lo[0]
            loy = ctile.lo[1]
            hix = ctile.hi[0]
            hiy = ctile.hi[1]
            if lox > vertx and hiy > verty:
                hixlm.append(lox)
            if loy > verty and hix > vertx:
                hiylm.append(loy)
        xhi = min(hixlm)
        yhi = min(hiylm)
        return [xhi, yhi]

    def get_enclosed_points(self, points, lo, hi):
        # Return list of points within [lo, hi]
        inpts = []
        lox = lo[0]
        loy = lo[1]
        hix = hi[0]
        hiy = hi[1]
        for pt in points:
            if (pt.x >= lox and pt.x <= hix and
                pt.y >= loy and pt.y <= hiy):
                inpts.append(pt)
        return inpts

    def get_dpt(self, a, b):
        return np.sqrt((a.x - b.x)**2 + (a.y - b.y)**2)

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
            
        # Start at lo
        # Find hi_prime giving minimum number of points
        # Expand 1 point at a time by points closest to lo in each direction until the next expansion gives geom_norm_resd > nres_threshold
        # Set hi by bisecting the distance between the most distal point such that geom_norm_resd <= nres_threshold in each direction
        # Create tile and repeat on the remaining space in the domain.
        # If domain is not covered, ERROR.
        
    def determine_subdivide(self):
        if len(self.points) == 3:
            return False
        elif len(self.points) < 3:
            print('ERROR: Tile contains less than 3 points')
            exit()
        p = Plane(self.points)
        if p.geom_norm_resd > self.nres_threshold:
            self.subdivide()

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

# Form upper level Tile to start
domain = Tile(pointlist, lo, hi)

# Recursively subdivide domain



