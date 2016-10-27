"""
Tile a 2-D domain containing Point objects depending on a decision function.

Donald Willcox
"""
import numpy as np
from Plane_nd import Plane_nd

class BCTypes(object):
    # Boundary Types
    up    = +1
    none  = None
    down  = -1

class DMCycle(object):
    # Dimension Cycler
    def __init__(self, dm=None):
        if not dm:
            print('ERROR: Must provide number of dimensions to DMCycle!')
            exit()
        self.dm = dm
        self.dims = range(self.dm)

    def cycle(self):
        t = self.dims.pop(0)
        self.dims.append(t)
        return t

class Point(object):
    def __init__(self, r=[], v=None):
        # Position in n-D space
        self.r = np.array(r)
        # Value of point (Scalar)
        self.v = v
        dm = len(r)
        # n-D Mask: Boundary Edge represented by Point
        self.bedge = [None for i in range(dm)]
        # n-D Mask: Boundary Type represented by Point
        self.btype = [BCTypes.none for i in range(dm)]

    def dist_to_pt(self, b):
        # Get distance between this Point and Point b
        dr = np.array(self.r) - np.array(b.r)
        return np.sqrt(np.sum(dr**2))

    def select_nn(self, plist):
        # Select the nearest neighbor from point list plist
        dmin = self.dist_to_pt(plist[0])
        p_nn = plist[0]
        pi_nn = 0
        for pi, p in enumerate(plist):
            dp = self.dist_to_pt(p)
            if dp < dmin:
                dmin = dp
                p_nn = p
                pi_nn = pi
        return pi_nn, p_nn


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
        self.points = points
        self.lo = lo
        self.hi = hi
        self.dm = dm
        self.nres_threshold = 0.5
        self.min_contain_points = self.dm + 1

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

class Domain(object):
    def __init__(self, points=[], lo=[], hi=[]):
        # The Domain is just a set of Point objects
        # and functions for tiling them into a set of Tile objects.        
        self.tiles = []
        self.lo = lo
        self.hi = hi
        self.dm = None
        self.points = points

        if lo and hi and len(lo) != len(hi):
            print('ERROR: lo and hi supplied with incongruous dimensions.')
            exit()
        else:
            self.dm = len(lo)
        
    def get_distal_point(self, refpt, points):
        # Get the most distant point from refpt among points
        dmax = 0.0
        pdst = None
        for p in points:
            dp = refpt.dist_to_pt(p)
            if dp > dmax:
                dmax = dp
                pdst = p
        return pdst

    def bc_init_mask_points(self):
        # Set initial boundary masks for points in domain
        # Initial because uses domain lo, hi
        dm_hi_pts = [None for i in range(self.dm)]
        dm_hi_val = self.lo
        dm_lo_pts = [None for i in range(self.dm)]
        dm_lo_val = self.hi
        for p in self.points:
            for di in range(self.dm):
                # Check high bc
                if p.r[di] > dm_hi_val[di]:
                    # Point is high, use point as bc
                    dm_hi_val[di] = p.r[di]
                    dm_hi_pts[di] = p
                if p.r[di] < dm_lo_val[di]:
                    # Point is low, use point as bc
                    dm_lo_val[di] = p.r[di]
                    dm_lo_pts[di] = p
        # Mask points in upper and lower bc lists
        for di, p in enumerate(dm_hi_pts):
            p.bedge[di] = self.hi[di]
            p.btype[di] = BCTypes.up
        for di, p in enumerate(dm_lo_pts):
            p.bedge[di] = self.lo[di]
            p.btype[di] = BCTypes.down

    def do_domain_tiling(self):
        dm_cycle = DMCycle(self.dm)
        # Get a dimension
        di = dm_cycle.cycle()
        # Select lower boundary point
        p_start = None
        pi_start = None
        for pi, p in enumerate(self.points):
            if p.btype[di] == BCTypes.down:
                p_start = p
                pi_start = pi
                break
        if not p_start and self.points:
            print('ERROR: Points remain but no boundary could be found!')
            exit()
        # Form tile with starting point
        t_start = Tile(points=[p_start], lo=p_start.r, hi=p_start.r, dm=self.dm)

        extend_cycle = DMCycle(self.dm)
        scratch_points = self.points[:]
        scratch_points.pop(pi)
        
        # Extend to n NN (not necessary to check planarity)
        # get a NN
        num_nn = 0
        while num_nn < self.dm+1:
            pi_nn, p_nn = p_start.select_nn(scratch_points)
            # Check to see if adding this NN will intersect existing Tiles.
            # DON

        # You have the n+1 point region R
        
        # extend R in each dimension
        # cut R out of domain (remove points and update point masks)
        # repeat while points remain in domain
    
# # Read Data
# dfname = 'output2.csv'
# raw_data = np.genfromtxt(dfname, delimiter=',', skip_header=1)
# # Each element of data is a row from the csv file, so convert to columns
# data = np.transpose(raw_data)
# # data[0] = Blocker factors
# xvec = data[0]
# # data[1] = Reimers factors
# yvec = data[1]
# # data[2] = CO WD Mass
# zvec = data[2]

# # Create list of Points
# pointlist = []
# for x, y, z in zip(xvec, yvec, zvec):
#     p = Point(x,y,z)
#     pointlist.append(p)

# # Get bounds on the x,y domain
# lo = [np.amin(xvec), np.amin(yvec)]
# hi = [np.amax(xvec), np.amax(yvec)]

# # Form Domain
# dom = Domain(pointlist, lo, hi)

        
