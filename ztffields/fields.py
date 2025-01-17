"""
fields
------

Module to handle ZTF Fields.

.. autosummary::

   Fields
   FieldVerts

"""
import pandas
import numpy as np
from shapely import geometry

from .utils import rot_xz_sph, _DEG2RA

# To access files
import os

_SOURCE_DATA = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
_FIELD_DATAFRAME = pandas.read_csv( os.path.join(_SOURCE_DATA, "ztf_fields.txt"), index_col="ID")


__all__ = ["Fields", "get_fieldid"]


def get_fieldid(grid=None, dec_range=None, ra_range=None, 
                gall_range=None, galb_range=None, 
                ecll_range=None, eclb_range=None, 
                ebv_range=None, verbose=False):
                
    """ get the list of fieldid parameters fullfilling the conditions

    Parameters
    ----------
    grid: Nonen, str
        Select the grid you want. If None the both will be considered.
        Available
        - 'main'
        - 'secondary'

    dec_range, ra_range, galb_range, gall_range, ecllg_range, eclb_range: None, list
        cut paramerers (in degree)
        3 formats available:
        - None: means no selection
        - [min,max]: means range to be considered, None means no limit. 
          example: decrange=[-10,None] means dec>-10.
        - [[min1,max1],[min2,max2], ...]: means zones to be considered.
          example: decrange=[[None, -10],[5,None]] will simply exclude the [-10,5] dec range.
        - None means no selection.

    ebv_range: None, list
        same format as e.g., 'dec_range' but for Mily way E(B-V) extinction.        
        
    Returns
    -------
    list
        list of int (fieldid)
    """
    def get_query(key, krange=None, merging_logic="or"):
        """ """
        def _build_2d_(kmin,kmax):
            """ """
            if kmin is None and kmax is None:
                return []
            if kmin is None:
                return [f"{key}<={kmax}"]
            elif kmax is None:
                return [f"{key}>={kmin}"]
            else:
                return [f"{kmin}<={key}<={kmax}"]
        
        if krange is None:
            return []
        if np.shape(krange) == (2,):
            return _build_2d_(*krange)
        if np.shape(krange) == (2,2):
            query = [_build_2d_(*krange_) for krange_ in krange]
            return ["("+f" {merging_logic} ".join(np.squeeze(query))+")"]
        
        raise ValueError(f"Cannot for the format of the input krange: {krange}")

    query = []
    # Grid
    gridid = None if (grid is None or grid in ["*","all"]) else get_grid_field(grid)
    query.append([] if gridid is None else ["ID in @gridid"])
    # Ra
    
    query.append(get_query("RA", ra_range))
    # Dec
    query.append(get_query("Dec", dec_range))
    # gall
    query.append(get_query("GalLong", gall_range))
    # galb
    query.append(get_query("GalLat", galb_range))
    # ecll
    query.append(get_query("EclLong", ecll_range))
    # eclb
    query.append(get_query("EclLat", eclb_range))
    # MW ebmv
    query.append(get_query("Ebv", ebv_range))
    
    queries = np.concatenate(query)
    if len(queries) == 0:
        fieldid = _FIELD_DATAFRAME.index
    else:
        fieldid = _FIELD_DATAFRAME.query( " & ".join(queries) ).index
        
    return list(fieldid)

def get_grid_field(which):
    """ """
    fieldid = _FIELD_DATAFRAME.index.values
    if which in ["main","Main","primary"]:
        return fieldid[fieldid<880]
    if which in ["aux","secondary", "auxiliary"]:
       return fieldid[fieldid>999]
    if which in ["all","*","both"]:
        return fieldid
        
    raise ValueError(f"Cannot parse which field grid you want {which}")


def get_layout(level):
    """ """
    if level.lower() == "ccd":
        layout = pandas.read_csv( os.path.join(_SOURCE_DATA, "ztf_ccd_layout.tbl")
                              ).rename(columns={"CCD ": "CCD"}) # strip
    elif level.lower() == "quadrant":
        layout = pandas.read_csv( os.path.join(_SOURCE_DATA, "ztf_ccd_quad_layout.tbl"))
        layout["EW"] *= -1 # error on the file
    else:
        raise NotImplementedError(f"only ccd and quadrant layout implemented. {level} given")
    
    return layout
        
        
class Fields( object ):
    """ Interact or access ZTF Fields. """
    
    FIELDS_RADEC = _FIELD_DATAFRAME[["RA","Dec"]].rename({"RA":"ra", "Dec":"dec"}, axis=1)
    CCD_COORDS = get_layout("ccd")
    QUAD_COORDS = get_layout("quadrant")


    def __init__(self, load_level=None):
        """ """
        if load_level is not None:
            self.load_level_geometry(load_level)

    def load_level_geometry(self, level):
        """ """
        level = np.atleast_1d(level)

        if "focalplane" in level or "*" in level:
            self._geodf_focalplane = self.get_field_geometry(level="focalplane")
        if "ccd" in level or "*" in level:
            self._geodf_ccd = self.get_field_geometry(level="ccd")
        if "quadrant" in level or "*" in level:
            self._geodf_quadrant = self.get_field_geometry(level="quadrant")

    def get_level_geodf(self, level):
        """ """
        if level=="focalplane":
            return self.geodf_focalplane
        
        if level=="ccd":
            return self.geodf_ccd
        
        if level=="quadrant":
            return self.geodf_quadrant
        
        raise ValueError(f"{level} not implemented.")
    
    # =============== #
    #   Methods       #
    # =============== #
    @classmethod
    def get_field_centroids(cls, fieldid=None):
        """ get field central coordinates [ra,dec]
        
        Parameters
        ----------
        fieldid: int, list
            name of the field (or list of).
            If None, no fieldid selection.
            
        Returns
        -------
        pandas.DataFrame
            dataframe (index as fieldid ; columns ra, dec)
        """
        radecs = cls.FIELDS_RADEC[["ra","dec"]]#.loc[fieldid]
        if fieldid is not None:
            radecs = radecs.loc[fieldid]
            
        return radecs
            
    @classmethod
    def get_field_geometry(cls, fieldid=None, level="ccd", **kwargs):
        """ get a geopandas dataframe for the fields.
        
        Parameters
        ----------
        fieldid: int, list
            name of the field (or list of).
            If None, no fieldid selection.

        level: str
            level of description of the camera.
            - focalplane: 1 polygon for the whole footprint
            - ccd: 1 polygon per CCD (16 per footprint then)
            - quadrant: 1 polygon per quadrant (64 per per footprint then)
            ccd and quadrant level will account for gaps between the CCDs.
            
        **kwargs goes to get_field_vertices
        
        Returns
        -------
        geopandas.GeoDataFrame
            geopandas GeoDataFrame containing the ['geometry'] columns.
            If level is ccd or quadrant, the dataframe has two level
            of indexes (fieldid, ccdid/rcid) ; only fieldid if level is focalplane.
            
        See also
        --------
        get_field_vertices: get the verticies (or polygon) for the given fields
        """
        from geopandas import geodataframe
        if fieldid is not None:
            fieldid = np.atleast_1d(fieldid) # needed for the following flatten()
            
        index, footprints = cls.get_field_vertices(fieldid=fieldid, 
                                                   level=level,
                                                   as_polygon=True)
        df = pandas.DataFrame(np.atleast_1d(footprints).flatten(), columns=["geometry"])
        df = df.join(index).set_index(index.columns.values.tolist())
        return geodataframe.GeoDataFrame(df)
            
    @classmethod
    def get_field_vertices(cls, fieldid=None, 
                            level="ccd",
                            inrad=False, steps=5, 
                            as_polygon=False):
        """ get the verticies (or polygon) for the given fields
        
        Parameters
        ----------
        fieldid: int, list
            name of the field (or list of).
            If None, no fieldid selection.

        level: str
            level of description of the camera.
            - focalplane: 1 polygon for the whole footprint
            - ccd: 1 polygon per CCD (16 per footprint then)
            - quadrant: 1 polygon per quadrant (64 per per footprint then)
            ccd and quadrant level will account for gaps between the CCDs.
            
        inrad: bool
            should the coordinates be given in radian ?
            False means in degree.
        
        steps: int
            number of points between two corners. These
            are used to account for square deformation.
            
        as_polygon: bool
            should this return a list of shapely.polygon rather than
            a list of (N,2) verticies ?
            
        Returns
        -------
        dataframe, list
            - the indexes (fieldid or (fieldid,ccdid/rcid) dependning on level)
            - the list of verticities ([steps, 2] or shapely.Polygon ; see as_polygon).
            
        See also
        --------
        get_field_geopandas: get a geopandas dataframe for the fields.
        """        
        radecs = cls.get_field_centroids(fieldid)
        ra, dec = radecs.values.T
        fieldids = radecs.index.values
        
        # focalplane level
        if level is None or level in ["focalplane","field"]:
            func = cls.get_focalplane_contours
            index = fieldids
            index = pandas.DataFrame(index, columns=["fieldid"])
        # ccd level
        elif level in ["ccd","ccdid"]:
            func = cls.get_ccd_contours
            index = [[fieldid,ccdid] for fieldid in fieldids for ccdid in np.arange(1,17)]
            index = pandas.DataFrame(index, columns=["fieldid","ccdid"])
        # rcid level            
        elif level in ["quadrant","rcid","quad"]:
            func = cls.get_quadrant_contours
            index = [[fieldid,rcid] for fieldid in fieldids for rcid in np.arange(0,64)]
            index = pandas.DataFrame(index, columns=["fieldid","rcid"])
        else:
            raise ValueError(f"level {level} is not available. use: focalplane, ccd, quadrant")
            
        footprints = func(ra=ra, dec=dec, inrad=inrad, steps=steps, as_polygon=as_polygon)
        return index, footprints

    
    @classmethod
    def get_field_multipolygon(cls, fieldid, level="ccd", **kwargs):
        """ get the a (multi)polygon for the given *single* field
        
        Parameters
        ----------
        fieldid: int
            name of the field.

        level: str
            level of description of the camera.
            - focalplane: 1 polygon for the whole footprint
            - ccd: 1 polygon per CCD (16 per footprint then)
            - quadrant: 1 polygon per quadrant (64 per per footprint then)
            ccd and quadrant level will account for gaps between the CCDs.

        **kwargs goes to get_field_geometry
        
        Returns
        -------
        shapely.(Multi)Polygon
            polygon if level is focalplane.            
        """
        # works only for single fieldid case
        fieldid = int(fieldid)
        footprint_df = cls.get_field_geometry(fieldid=fieldid, level=level, **kwargs)
        if level == "focalplane":
            return footprint_df.loc[fieldid].values[0]
        
        poly_list = footprint_df.xs(fieldid)["geometry"].values.tolist()

        return geometry.MultiPolygon(poly_list)
    
    # ---------- #
    #  Corners   #
    # ---------- #
    @classmethod
    def get_contours(cls, level="quadrant", inrad=False, steps=5,
                         ra=0, dec=0, as_polygon=False,
                         allow_multipolygon=True, **kwargs):
        """ """
        func = getattr(cls,f"get_{level}_contours")
        contours = func(inrad=inrad, steps=steps, ra=ra, dec=dec, as_polygon=as_polygon, **kwargs)
        if as_polygon and level in ["ccd", "quadrant"] and allow_multipolygon:
            contours = geometry.MultiPolygon(list(contours))

        return contours
    
    # Quadrant    
    @classmethod    
    def get_quadrant_contours(cls, rcid=None, 
                            inrad=False, steps=5, 
                            ra=0, dec=0, as_polygon=False):
        """ get the quadrant contours.
        
        Parameters
        ----------
        rcid: int, list
            id of the quadrant (or list of).
            If None, all 64 quadrants will be used.
            
        inrad: bool
            should the coordinates be given in radian ?
            False means in degree.
        
        steps: int
            number of points between two corners. These
            are used to account for square deformation.
            
        ra: float, array
            field ra centroid (in degree if not inrad)
            
        dec: float, array
            field dec centroid
            
        as_polygon: bool
            should the return a shapely polygon or the vertices (N,2) array
            
        Returns
        -------
        array or lit of polygon
            see as_polygon
        
        See also
        --------
        get_ccd_contours: external countour of a full ccd
        get_focalplane_contours: external countour of the whole focal plane (Field of view)
        """
        upper_left_corner = cls.QUAD_COORDS.groupby("Quad").max()
        lower_right_corner = cls.QUAD_COORDS.groupby("Quad").min()
        if rcid is not None:
            upper_left_corner = upper_left_corner.loc[rcid]
            lower_right_corner = lower_right_corner.loc[rcid]
            
        return cls._derive_corners_(upper_left_corner, lower_right_corner, 
                                    ra=ra, dec=dec,
                                    inrad=inrad, steps=steps,
                                    as_polygon=as_polygon)
    # CCD
    @classmethod    
    def get_ccd_contours(cls, ccdid=None, 
                            inrad=False, steps=5, 
                            ra=0, dec=0, as_polygon=False):
        """ get the full ccd contours.
        
        Parameters
        ----------
        ccdid: int, list
            id of the ccd (or list of).
            If None, all 16 quadrants will be used.
            
        inrad: bool
            should the coordinates be given in radian ?
            False means in degree.
        
        steps: int
            number of points between two corners. These
            are used to account for square deformation.
            
        ra: float, array
            field ra centroid (in degree if not inrad)
            
        dec: float, array
            field dec centroid
            
        as_polygon: bool
            should the return a shapely polygon or the vertices (N,2) array
            
        Returns
        -------
        array or lit of polygon
            see as_polygon
        
        See also
        --------
        get_quadrant_contours:  countour of a single quadrant (ccd amplifier)
        get_focalplane_contours: external countour of the whole focal plane (Field of view)
        """
        upper_left_corner = cls.CCD_COORDS.groupby("CCD").max()
        lower_right_corner = cls.CCD_COORDS.groupby("CCD").min()
        if ccdid is not None:
            upper_left_corner = upper_left_corner.loc[ccdid]
            lower_right_corner = lower_right_corner.loc[ccdid]
            
        return cls._derive_corners_(upper_left_corner, lower_right_corner, 
                                    ra=ra, dec=dec,
                                    inrad=inrad, steps=steps,
                                    as_polygon=as_polygon)
    
    # Focal Plane
    @classmethod
    def get_focalplane_contours(cls, inrad=False, steps=5, 
                               ra=0, dec=0,
                               as_polygon=False):
        """ get the whole focal plane contours (field of view)
        
        Parameters
        ----------            
        inrad: bool
            should the coordinates be given in radian ?
            False means in degree.
        
        steps: int
            number of points between two corners. These
            are used to account for square deformation.
            
        ra: float, array
            field ra centroid (in degree if not inrad)
            
        dec: float, array
            field dec centroid
            
        as_polygon: bool
            should the return a shapely polygon or the vertices (N,2) array
            
        Returns
        -------
        array or lit of polygon
            see as_polygon
        
        See also
        --------
        get_quadrant_contours:  countour of a single quadrant (ccd amplifier)
        get_ccd_contours: external countour of a full ccd
        """
        upper_left_corner = cls.CCD_COORDS.max()
        lower_right_corner = cls.CCD_COORDS.min()
    
        return cls._derive_corners_(upper_left_corner, lower_right_corner, 
                                    ra=ra, dec=dec,
                                    inrad=inrad, steps=steps,
                                    as_polygon=as_polygon)
    
    # =============== #
    #    Internal     #
    # =============== #
    

    @staticmethod
    def _verts_to_polygon_(verts):
        """ """
        shapes = np.shape(verts)
        if len(shapes)==2:
            fields_countours = geometry.Polygon(verts)
        else:
            shape_flat = tuple([np.prod(shapes[:-2])] + list(shapes[-2:]))
            fields_countours = np.reshape([geometry.Polygon(f_) 
                                           for f_ in verts.reshape(shape_flat)], shapes[:-2])
        return fields_countours
    
    @classmethod
    def _derive_corners_(cls, 
                         upper_left_corner, lower_right_corner, 
                         ra, dec, steps=5, squeeze=True, 
                         inrad=False, as_polygon=False):
        """ """
        
        ewmin = -np.atleast_1d(upper_left_corner["EW"])
        nsmax = +np.atleast_1d(upper_left_corner["NS"])
        ewmax = -np.atleast_1d(lower_right_corner["EW"])
        nsmin = +np.atleast_1d(lower_right_corner["NS"])

        ra1  = (np.linspace(ewmax, ewmin, steps)/np.cos(nsmax*_DEG2RA)).T
        dec1 = (np.ones((steps,1))*nsmax).T
        #
        dec2  = np.linspace(nsmax,nsmin, steps).T        
        ra2   = ewmin[:,None]/np.cos(dec2*_DEG2RA)
        #
        ra3 = (np.linspace(ewmin,ewmax, steps)/np.cos(nsmin*_DEG2RA)).T
        dec3 = (np.ones((steps,1))*nsmin).T
        #
        dec4  = np.linspace(nsmin,nsmax, steps).T
        ra4 = ewmax[:,None]/np.cos(dec4*_DEG2RA)
        
        ra_bd = np.concatenate((ra1, ra2, ra3, ra4  ), axis=1)  
        dec_bd = np.concatenate((dec1, dec2, dec3, dec4 ), axis=1)

        ra_, dec_ = rot_xz_sph(np.moveaxis(ra_bd,0,1),
                               np.moveaxis(dec_bd,0,1),
                               np.moveaxis(np.atleast_3d(dec),0,1))
        ra_ += np.moveaxis(np.atleast_3d(ra), 0, 1)

        if inrad:
            ra_ *= _DEG2RA
            dec_ *= _DEG2RA

        radec = np.moveaxis([ra_,dec_], (0,1,2,3), (3,0,2,1))
        if squeeze:
            radec = np.squeeze(radec)

        if as_polygon:
            return cls._verts_to_polygon_(radec)
        return radec
    
    # =============== #
    #   Properties    #
    # =============== #
    @property
    def geodf_focalplane(self):
        """ field geometry at the focal plane level """
        if not hasattr(self,"_geodf_focalplane"):
            self.load_level_geometry("focalplane")
            
        return self._geodf_focalplane
    
    @property
    def geodf_ccd(self):
        """ field geometry at the ccd level """
        if not hasattr(self,"_geodf_ccd"):
            self.load_level_geometry("ccd")
            
        return self._geodf_ccd

    @property
    def geodf_quadrant(self):
        """ field geometry at the quadrant level """
        if not hasattr(self,"_geodf_quadrant"):
            self.load_level_geometry("quadrant")
            
        return self._geodf_quadrant


# ================ #
#  Fast Vertices   #
# ================ #
class FieldVerts( object ):
    
    def __init__(self, load_which="*"):
        """ """
        if load_which in ["all", "*"]:
            load_which = ["quadrant", "ccd", "focalplane"]
        self._fieldids = {}
        self._verts = {}
        
        # Focal Plane
        for which in np.atleast_1d(load_which):
            fieldid, verts = Fields.get_field_vertices(level=which)
            self._fieldids[which] = fieldid
            self._verts[which] = verts

    # ============= #
    #  Methods      #
    # ============= #
    def get_centroid(self, level, as_dataframe=False):
        """ """
        verts = self.verts.get(level)
        if verts is None:
            raise ValueError(f"no vertices for {level=}")

        verts_npoints = verts.shape[-2] # like 20
        # as numpy
        centroids = np.mean(verts.reshape(np.prod(verts.shape[:-2]),
                                          verts_npoints, 2), axis=1)
        
        if as_dataframe: # as dataframe
            centroids = pandas.DataFrame(centroids, columns=["ra", "dec"])
            if level == "focalplane":
                centroids = centroids.set_index(self.fieldids.get(level)["fieldid"])
            else:
                index = pandas.MultiIndex.from_frame(self.fieldids.get(level))
                centroids = centroids.set_index(index)
            
        return centroids

    def get_corner(self, level, which):
        """ 
        Parameters
        ----------
        level: str
            focalplane, ccd, or quadrant

        which: str, int
            which corner. 
            - 0: upper right
            - 1: upper left
            - 2: lower left
            - 3: lower right
        """
        # input parsing
        WHICH_DICT = {"upper right":0, "upper left":1,
                      "lower left":2, "lower right":3}
        if type(which) is str:
            which = WHICH_DICT.get(which) # now an int
            
        # grabbing vertices
        verts = self.verts.get(level)
        if verts is None:
            raise ValueError(f"no vertices for {level=}")

        # getting the correct vert entry given the number of steps        
        verts_nsteps = int(verts.shape[-2]/4) # like 20/4 = 5
        if level == "focalplane":
            corner =  verts[:, int(which*verts_nsteps), :]
        else:
            corner =  verts[:, :, int(which*verts_nsteps), :].reshape(
                np.prod(verts.shape[:-2]), 2)

        corner = pandas.DataFrame(corner, columns=["ra", "dec"])
        if level == "focalplane":
            corner = corner.set_index(self.fieldids.get(level)["fieldid"])
        else:
            index = pandas.MultiIndex.from_frame(self.fieldids.get(level))
            corner = corner.set_index(index)
            
        return corner
    
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def fieldids(self):
        """ """
        return self._fieldids

    @property
    def verts(self):
        """ """
        return self._verts
