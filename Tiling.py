"""
Tile an n-D domain containing Point objects depending on a decision function.

Donald Willcox
"""
import numpy as np
from Plane_nd import Plane_nd

class BCTypes(object):
    # Boundary Types
    up    = +1
    none  = None
    down  = -1
    tile  = +2
    point = +3

class DMCycle(object):
    # Dimension Cycler
    def __init__(self, dm=None):
        if not dm:
            print('ERROR: Must provide number of dimensions to DMCycle!')
            exit()
        self.dm = dm
        self.dims = [i for i in range(self.dm)]

    def cycle(self):
        t = self.dims.pop(0)
        self.dims.append(t)
        return t

class TilingError(Exception):
    """Error class for various kinds of tiling errors that can arise."""
    def __init__(self, err_tile=None, scratch_points=None, message=''):
        # err_tile is the Tile object we were attempting to extend
        self.err_tile = err_tile
        self.scratch_points = scratch_points
        self.message = message

class Point(object):
    def __init__(self, r=[], v=None):
        # Position in n-D space
        self.r = np.array(r)
        # Value of point (Scalar)
        self.v = v
        self.dm = len(r)
        # n-D Mask: Boundary Edge represented by Point
        self.bedge = [None for i in range(self.dm)]
        # n-D Mask: Boundary Type represented by Point
        self.btype = [BCTypes.none for i in range(self.dm)]

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
        
class Tile(object):
    def __init__(self, points=[], lo=[], hi=[], dm=None):
        self.points = points
        self.lo = lo
        self.hi = hi
        self.dm = dm
        if points and not dm:
            self.dm = points[0].dm
        if not lo and not hi:
            self.boundary_minimize()

    def get_volume(self):
        """
        Computes volume of the Tile. If [lo, hi] is undefined, return None.
        """
        if not self.lo or not self.hi:
            return None
        dr = np.array(self.hi) - np.array(self.lo)
        return np.prod(dr)

    def boundary_minimize(self):
        """
        Given the points in the Tile, set the boundary
        defined by [lo, hi] to the minimum surface enclosing the points.
        """
        if self.points:
            self.lo = self.points[0].r
            self.hi = self.points[0].r
            for p in self.points:
                self.lo = np.minimum(self.lo, p.r)
                self.hi = np.maximum(self.hi, p.r)

    def extend_points(self, plist=[]):
        """
        Given the list of points (plist), extends the Tile if
        necessary to enclose them.

        Set the Tile boundaries to the minimum volume enclosing
        the provided points.

        Do nothing if no points are provided.
        """
        if not plist:
            return
        self.points += plist
        self.boundary_minimize()

    def overlaps_point_dimension(self, refpoint, di):
        """
        Checks to see if self overlaps refpoint in the dimension di:

        refpoint must be a Point object

        di must be an integer in range(self.dm)
        """
        refp_olap = True
        if (refpoint.r[di] >= self.hi[di] or
            refpoint.r[di] <= self.lo[di]):
            # tiles do not overlap
            refp_olap = False
        return refp_olap

    def overlaps_tile_dimension(self, reftile, di):
        """
        Checks to see if self overlaps reftile in the dimension di:

        reftile must be a Tile object

        di must be an integer in range(self.dm)
        """
        reft_olap = True
        if (reftile.lo[di] >= self.hi[di] or
            reftile.hi[di] <= self.lo[di]):
            # tiles do not overlap
            reft_olap = False
        return reft_olap

    def overlaps_tiles(self, tlist=[]):
        """
        Checks to see if self overlaps any of the tiles in tlist

        Returns list of tiles in tlist which overlap self

        Returns the empty list if no tiles in tlist overlap self
        """
        if not tlist:
            return False
        olap = []
        for reft in tlist:
            # Check overlap between self and reft
            reft_olap = True
            for di in range(self.dm):
                reft_olap = self.overlaps_tile_dimension(reft, di)
                if not reft_olap:
                    # tiles do not overlap
                    break
            if reft_olap:
                olap.append(reft)
        return olap        

    def extend_min_volume(self, plist=[], avoid_tiles=None, decision_fun=None):
        """
        Given the list of points (plist), extends the Tile
        by adding one point from plist to Tile where the point 
        is selected from plist such that it minimizes the volume
        of Tile.

        Returns plist where the selected point is popped from the list.

        If avoid_tiles is passed, it should be a list of Tile objects.
        The current Tile will then only be extended such that it does
        not intersect the tiles in avoid_tiles.

        If a function is passed as decision_fun, this Tile will be passed to
        the 'decision function' to determine whether to extend the tile.
        decision_fun should take a single Tile argument and return True or False
        """
        min_vol_point_i = None
        min_vol_point = None
        min_vol = None
        for i, p in enumerate(plist):
            pext = self.points + [p]
            stile = Tile(points=pext)
            svol = stile.get_volume()
            dbool = True
            if callable(decision_fun):
                dbool = decision_fun(stile)
            if (((not min_vol_point) or svol < min_vol)
                and
                dbool
                and
                not stile.overlaps(avoid_tiles)):
                min_vol = svol
                min_vol_point = p
                min_vol_point_i = i
        if not min_vol_point:
            # If the above could find no point to extend, then do nothing
            # Return point list and False, indicating no extension
            return plist, False
        else:
            # Else, extend this Tile
            self.extend_points([min_vol_point])
            # Return reduced point list and True, indicating extension
            return plist.pop(i), True

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

        # Set up dimension cycling
        self.dm_cycle = DMCycle(self.dm)
        
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

    def set_tile_boundaries(self, atile):
        """
        Given atile, sets its [lo, hi] boundaries in each dimension.

        Also updates the boundary masks for adjacent points
        in the tiling list scratch_points.
        """
        # Expand Tile in each dimension as possible
        for di in range(self.dm):
            # Find the tiles which atile does not overlap in dimension di
            # but does overlap in every other dimension.
            # These tiles set the bounds on extensions along dimension di.
            # (If there is an additional dimension along which the tiles
            # do not overlap, then no constraint can be made along di.)
            otiles = []
            for ktile in self.tiles:
                if atile.overlaps_tile_dimension(ktile, di):
                    # ktile overlaps along di, so can't constrain di
                    continue
                kandidate = True
                for dj in range(self.dm):
                    if dj==di:
                        continue
                    if not atile.overlaps_tile_dimension(ktile, dj):
                        # ktile doesn't overlap along dj, dj =/= di
                        # so can't constrain di
                        kandidate = False
                        break
                if kandidate:
                    otiles.append(ktile)

            # Find the points which atile does not overlap in dimension di
            # but does overlap in every other dimension.
            # These points set the bounds on extensions along dimension di.
            # (If there is an additional dimension along which tile and points
            # do not overlap, then no constraint can be made along di.)
            opoints = []
            for p in self.points:
                if atile.overlaps_point_dimension(p, di):
                    # p overlaps along di, so can't constrain di
                    continue
                kandidate = True
                for dj in range(self.dm):
                    if dj==di:
                        continue
                    if not atile.overlaps_point_dimension(p, dj):
                        # p doesn't overlap along dj, dj =/= di
                        # so can't constrain di
                        kandidate = False
                        break
                if kandidate:
                    opoints.append(p)
            
            # Setup bc data structures for figuring out boundaries
            lo_bc = BCTypes.none
            lo_bc_type = BCTypes.none
            lo_bc_object = BCTypes.none
            hi_bc = BCTypes.none
            hi_bc_type = BCTypes.none
            hi_bc_object = BCTypes.none

            # Get tile constraint on di for [lo, hi]
            for btile in otiles:
                # Check if btile can constrain lo along di
                if btile.hi[di] <= atile.lo[di]:
                    if lo_bc == BCTypes.none or btile.hi[di] > lo_bc:
                        lo_bc = btile.hi[di]
                        lo_bc_type = BCTypes.tile
                        lo_bc_object = btile

                # Check if btile can constrain hi along di
                if btile.lo[di] >= atile.hi[di]:
                    if hi_bc == BCTypes.none or btile.lo[di] < hi_bc:
                        hi_bc = btile.lo[di]
                        hi_bc_type = BCTypes.tile
                        hi_bc_object = btile

            # Get point constraint on di for [lo, hi]
            for p in opoints:
                # Check if p can constrain lo along di
                if p.r[di] <= atile.lo[di]:
                    phalf = 0.5*(p.r[di] + atile.lo[di])
                    if lo_bc == BCTypes.none or phalf > lo_bc:
                        lo_bc = phalf
                        lo_bc_type = BCTypes.point
                        lo_bc_object = p

                # Check if p can constrain hi along di
                if p.r[di] >= atile.hi[di]:
                    phalf = 0.5*(p.r[di] + atile.hi[di])
                    if hi_bc == BCTypes.none or phalf < hi_bc:
                        hi_bc = phalf
                        hi_bc_type = BCTypes.point
                        hi_bc_object = p
            
            # If neither point nor tile constraint, use domain [lo, hi]
            if lo_bc == BCTypes.none:
                lo_bc = self.lo[di]
            if hi_bc == BCTypes.none:
                hi_bc = self.hi[di]

            # If point constraint, make that point a boundary along di
            if lo_bc_type == BCTypes.point:
                # Set point bc masks wrt points remaining in domain
                lo_bc_object.bedge[di] = lo_bc
                lo_bc_object.btype[di] = BCTypes.upper
            if hi_bc_type == BCTypes.point:
                # Set point bc masks wrt points remaining in domain
                hi_bc_object.bedge[di] = hi_bc
                hi_bc_object.btype[di] = BCTypes.lower
            # Go to next dimension

    def form_tile(self, gnr_thresh=0.5):
        # Cycle through dimensions
        di = self.dm_cycle.cycle()
        
        # Select lower boundary point in dimension di
        p_start = None
        pi_start = None
        for pi, p in enumerate(self.points):
            if p.btype[di] == BCTypes.down or p.btype[di] == BCTypes.up:
                p_start = p
                pi_start = pi
                break
        if not p_start and self.points:
            print('ERROR: Points remain but no boundary could be found!')
            print('Dimension {}'.format(di))
            print('Number of Domain Points {}'.format(len(self.points)))
            print('Number of Domain Tiles {}'.format(len(self.tiles)))
            exit()
            
        # Form tile with starting point
        atile= Tile(points=[p_start], lo=p_start.r, hi=p_start.r, dm=self.dm)
        self.scratch_points.pop(pi_start)
        
        # Extend to enclose a total n+1 points for n-D parameter space.
        # More points may be enclosed if exactly n+1 isn't possible
        canex = True
        while len(atile.points) < self.dm+1 and canex:
            self.scratch_points, canex = atile.extend_min_volume(plist=self.scratch_points,
                                                                 avoid_tiles=self.tiles)
            
        # Check the number of points, if it's less than n+1, raise an exception!
        if len(atile.points) < self.dm+1:
            raise TilingError(atile, self.scratch_points, 
                              'Could not enclose n+1 points!')
            
        # extend Tile checking the fit decision function dfun
        dfun  = (lambda t: t.get_geom_norm_resd() < gnr_thresh)
        canex = True
        while canex:
            self.scratch_points, canex = atile.extend_min_volume(plist=self.scratch_points,
                                                                 avoid_tiles=self.tiles,
                                                                 decision_fun=dfun)
            
        # Check that at least n+1 points remain in the domain, otherwise raise an exception!
        if len(self.scratch_points) < self.dm+1:
            raise TilingError(atile, self.scratch_points, 'Fewer than n+1 points remain!')
            
        # set boundaries of atile and update point boundary masks
        self.set_tile_boundaries(atile)

        # Add atile to tiles in this domain
        self.tiles.append(atile)

    def do_domain_tiling(self, gnr_thresh=0.5):
        # Initialize a list of scratch points for tiling
        self.scratch_points = self.points[:]
        # Clear current tiling
        self.tiles = []
        # Tile the domain given gnr_thresh
        while self.scratch_points:
            try:
                self.form_tile(gnr_thresh)
            except TilingError as terr:
                print(terr.message)
                print('Number of points in attempted tile: {}'.format(
                    len(terr.err_tile.points)))
                print('Number of points remaining in domain: {}'.format(
                    len(terr.scratch_points)))
                raise
