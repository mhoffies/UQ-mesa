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

class TETypes(object):
    """Tiling Error Types"""
    # Tiling cannot enclose enough points to constrain the fit
    cannot_enclose_enough_points = 1
    # Too few points remain in domain to constrain a fit on a new tile
    few_points_remain = 2

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
    def __init__(self, err_type, err_tile=None, scratch_points=None, message=''):
        # err_tile is the Tile object we were attempting to extend
        self.err_type = err_type
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
    def __init__(self, points=None, fit_guess=None, dm=None):
        # fit_guess should provide, well, a guess for the fit parameters
        # if fit_guess isn't supplied, it will be estimated from the points.
        self.cpars = None # Length n+1 for n-D space
        self.dm = dm
        self.resd  = None
        self.norm_resd = None
        self.geom_norm_resd = None
        self.compute_pars(points, fit_guess)
        
    def compute_pars(self, points, fit_guess):
        if not points:
            return
        ivars = np.array([p.r for p in points])
        dvars = np.array([p.v for p in points])
        fitter = Plane_nd(ivars, dvars, self.dm)
        popt, pcov = fitter.dolsq(fit_guess)
        dpfit = np.array([fitter.fplane(ivr, popt) for ivr in ivars])
        self.resd = dvars - dpfit
        self.norm_resd = self.resd/dvars
        self.geom_norm_resd = np.sqrt(np.sum(self.norm_resd**2))
        
class Tile(object):
    def __init__(self, points=[], lo=[], hi=[], fit_guess=None, dm=None):
        self.points = points
        self.fit_guess = fit_guess
        self.lo = lo
        self.hi = hi
        self.dm = dm
        self.geom_norm_resd = None
        if all(points) and not dm:
            self.dm = points[0].dm
        if not list(self.lo) or not list(self.hi):
            self.boundary_minimize()

    def print_tile_report(self):
        """Prints report of this Tile"""
        print('---TILE REPORT---')
        print('---GEOM. MEAN NORM RESD. = {} ---'.format(self.get_geom_norm_resd()))
        print('-------POINTS------')
        for p in self.points:
            print('{}: {}'.format(p.r, p.v))

    def get_volume(self):
        """
        Computes volume of the Tile. If [lo, hi] is undefined, return None.
        """
        if not list(self.lo) or not list(self.hi):
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
            stile = Tile(points=pext, fit_guess=self.fit_guess)
            svol = stile.get_volume()
            dbool = True
            if callable(decision_fun):
                dbool = decision_fun(stile)
            if (((not min_vol) or svol < min_vol)
                and
                dbool
                and
                not stile.overlaps_tiles(avoid_tiles)):
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
            plist.pop(min_vol_point_i)
            return plist, True

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
        p = Plane(points=self.points, fit_guess=self.fit_guess,
                  dm=self.dm)
        self.fit_guess = p.cpars
        self.geom_norm_resd = p.geom_norm_resd
        return self.geom_norm_resd

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

        # Set up boundary masks for points
        if all(points) and all(lo) and all(hi):
            self.bc_init_mask_points(self.points)

    def print_domain_report(self):
        # Prints full (scratch) domain data
        print('---DOMAIN REPORT---')
        print('DOMAIN LO = {}'.format(self.lo))
        print('DOMAIN HI = {}'.format(self.hi))
        print('-------POINTS------')
        for p in self.scratch_points:
            print('{}: {}'.format(p.r, p.v))
        print('--------TILES------')
        for t in self.tiles:
            print('lo = {}, hi = {}, npts = {}, geom_norm_resd = {}'.format(t.lo, t.hi,
                                                                            len(t.points),
                                                                            t.get_geom_norm_resd()))
            t.print_tile_report()
        
    def bc_init_mask_points(self, plist):
        # Set initial boundary masks for list of points plist given domain boundaries
        # Initial because uses domain lo, hi
        dm_hi_pts = [None for i in range(self.dm)]
        dm_hi_val = self.lo[:]
        dm_lo_pts = [None for i in range(self.dm)]
        dm_lo_val = self.hi[:]
        for p in plist:
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
            print('BC HI PT {}'.format(di))
            print('position: {}'.format(p.r))
            p.bedge[di] = self.hi[di]
            p.btype[di] = BCTypes.up
            print('bedge: {}'.format(p.bedge[di]))
            print('btype: {}'.format(p.btype[di]))
            
        for di, p in enumerate(dm_lo_pts):
            print('BC LO PT {}'.format(di))
            print('position: {}'.format(p.r))
            p.bedge[di] = self.lo[di]
            p.btype[di] = BCTypes.down
            print('bedge: {}'.format(p.bedge[di]))
            print('btype: {}'.format(p.btype[di]))

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
            for p in self.scratch_points:
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
                        lo_bc = btile.hi[di] + np.finfo(float).eps
                        lo_bc_type = BCTypes.tile
                        lo_bc_object = btile

                # Check if btile can constrain hi along di
                if btile.lo[di] >= atile.hi[di]:
                    if hi_bc == BCTypes.none or btile.lo[di] < hi_bc:
                        hi_bc = btile.lo[di] - np.finfo(float).eps
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

            # Now implement [lo_bc, hi_bc] for this tile and dimension di
            atile.lo[di] = lo_bc
            atile.hi[di] = hi_bc

            # If point constraint, make that point a boundary along di
            if lo_bc_type == BCTypes.point:
                # Set point bc masks wrt points remaining in domain
                lo_bc_object.bedge[di] = lo_bc
                lo_bc_object.btype[di] = BCTypes.up
            if hi_bc_type == BCTypes.point:
                # Set point bc masks wrt points remaining in domain
                hi_bc_object.bedge[di] = hi_bc
                hi_bc_object.btype[di] = BCTypes.down

            # Explain Yourself!
            print('Setting Tile Boundaries along dimension {} for Reasons:'.format(di))
            print('lo reason: {}'.format(lo_bc_type))
            print('hi reason: {}'.format(hi_bc_type))
            # Go to next dimension

    def form_tile(self, gnr_thresh=None):
        # Cycle through dimensions
        di = self.dm_cycle.cycle()

        print('Executing Domain.form_tile()')
        print('Number of points in domain: {}'.format(len(self.scratch_points)))
        print('Number of tiles in domain: {}'.format(len(self.tiles)))

        # Select lower boundary point in dimension di
        p_start = None
        pi_start = None
        for pi, p in enumerate(self.scratch_points):
            if p.btype[di] == BCTypes.down or p.btype[di] == BCTypes.up:
                p_start = p
                pi_start = pi
                break
        if not p_start and self.scratch_points:
            print('Points remain but no boundary could be found!')
            print('Dimension {}'.format(di))
            print('Number of Domain Points {}'.format(len(self.scratch_points)))
            print('Number of Domain Tiles {}'.format(len(self.tiles)))
            print('Re-establishing boundaries')
            self.bc_init_mask_points(self.scratch_points)
            self.form_tile(gnr_thresh)
        else:
            print('Found starting boundary point')
            print('Start index: {}'.format(pi_start))
            print('Start position: {}'.format(p_start.r))
            
        # Form tile with starting point
        atile = Tile(points=[p_start], lo=p_start.r, hi=p_start.r, dm=self.dm)
        self.scratch_points.pop(pi_start)

        print('Getting at least n+1 points')
        # Extend to enclose a total n+1 points for n-D parameter space.
        # More points may be enclosed if exactly n+1 isn't possible
        canex = True
        while len(atile.points) < self.dm+1 and canex:
            self.scratch_points, canex = atile.extend_min_volume(plist=self.scratch_points,
                                                                 avoid_tiles=self.tiles)
            print('Attempted tile has {} points'.format(len(atile.points)))
        print('Obtained {} points'.format(len(atile.points)))
        
        # Check the number of points, if it's less than n+1, raise an exception!
        if len(atile.points) < self.dm+1:
            raise TilingError(err_type=TETypes.cannot_enclose_enough_points,
                              err_tile=atile,
                              scratch_points=self.scratch_points, 
                              message='Could not enclose n+1 points!')
            
        # extend Tile checking the fit decision function dfun
        print('Extending initial tile')
        if gnr_thresh:
            dfun = (lambda t: t.get_geom_norm_resd() < gnr_thresh)
        else:
            dfun = (lambda t: True)
        canex = True
        while canex:
            self.scratch_points, canex = atile.extend_min_volume(plist=self.scratch_points,
                                                                 avoid_tiles=self.tiles,
                                                                 decision_fun=dfun)
            print('Attempted tile has {} points'.format(len(atile.points)))
                        
        # set boundaries of atile and update point boundary masks
        print('Updating point boundary masks')
        self.set_tile_boundaries(atile)

        # Add atile to tiles in this domain
        print('Adding tile to domain')
        self.tiles.append(atile)

        # Check that at least n+1 points remain in the domain, otherwise warn user!
        if self.scratch_points and len(self.scratch_points) < self.dm+1:
            raise TilingError(err_type=TETypes.few_points_remain,
                              err_tile=atile,
                              scratch_points=self.scratch_points,
                              message='Fewer than n+1 points remain!')

    def extend_existing_tiles(self, gnr_thresh=None):
        """
        Extends all existing tiles, gobbling up scratch_points as possible.
        """
        canex_tiles = False
        for i, atile in enumerate(self.tiles):
            other_tiles = self.tiles[:]
            other_tiles.pop(i)
            
            if gnr_thresh:
                dfun  = (lambda t: t.get_geom_norm_resd() < gnr_thresh)
            else:
                dfun = (lambda t: True)
                
            print('Existing tile {} has {} points'.format(i, len(atile.points)))
            self.scratch_points, canex = atile.extend_min_volume(plist=self.scratch_points,
                                                                 avoid_tiles=other_tiles,
                                                                 decision_fun=dfun)

            canex_tiles = canex_tiles or canex
            # set boundaries of atile and update point boundary masks
            print('Updating point boundary masks')
            self.set_tile_boundaries(atile)
        return canex_tiles

    def do_domain_tiling(self, gnr_thresh=None):
        # Initialize a list of scratch points for tiling
        self.scratch_points = self.points[:]
        # Clear current tiling
        self.tiles = []
        # Tile the domain given gnr_thresh
        try:
            while self.scratch_points:
                self.form_tile(gnr_thresh)
            self.print_domain_report()
        except TilingError as terr:
            print(terr.message)
            print('Number of points in attempted tile: {}'.format(
                len(terr.err_tile.points)))
            print('Number of points remaining in domain: {}'.format(
                len(terr.scratch_points)))
            if (terr.err_type == TETypes.few_points_remain or
                terr.err_type == TETypes.cannot_enclose_enough_points):
                # Distribute remaining points among existing tiles
                # gnr_thresh constraint is ignored, so
                # tile overlap is the only constraint.
                canex = True
                while self.scratch_points and canex:
                    canex = self.extend_existing_tiles()
                self.print_domain_report()
            else:
                raise
