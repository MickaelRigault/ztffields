"""
projection
----------

Project targets (ra, dec) into the fields.

.. autosummary::

   radec_to_fieldid
   fieldid_to_radec
   get_camera_location
   spatialjoin_radec_to_fields
   FieldProjection

"""

import warnings
import pandas
import numpy as np
#
import geopandas
from shapely import geometry

#
from .fields import Fields, FieldVerts

__all__ = ["radec_to_fieldid", "fieldid_to_radec",
               "get_camera_location",
               "FieldProjection"]

FIELD_VERTICES = FieldVerts(load_which="all")
RA_SIZE, DEC_SIZE = 3080, 3072
# ============== #
#                #
#  shortcuts     #
#                #
# ============== #

def radec_to_fieldid(radecs, level="focalplane", **kwargs):
    """ """
    return FieldProjection.radec_to_fieldid(radecs, level=level, **kwargs)

def fieldid_to_radec(fieldid=None, level="focalplane", **kwargs):
    """ 
    Parameters
    ----------
    fieldid: int, list, array
        If None, all is returned.
        expected shape of fieldid depends on requested level:
        - focalplane: 1d array (N,) will return N entries e.g.
            (fieldid_1, fieldid_2, ...)
        - ccd, quadrant: 1 (N,) or 2d array (N, 2) will return N entries 
            if (N,): expected: (fieldid_1, 
                                fieldid_2, ...)
            else: expected: ( [fieldid_1,ccdid_1], [fieldid_1,ccdid_2], 
                              [fieldid_2,ccdid_3], ...)
                              
        # Caution example:
        fieldid_to_radec([427,14], level="ccd") 
        => returns 32 entries (all 16 CCDs for fields 427 and 14
        fieldid_to_radec([ [427,14] ], level="ccd") 
        => returns 1 entry, the CCD 14 of field 427
        
        # Alternative
        for more option, simply use fieldid = None and parse the 
        pandas multi-index dataframe after. 
        For instance, to get all 13 of all fields:
        fieldid_to_radec(level="ccd").xs(13, level=1).
        
    level: str
        level of description of the camera.
        - focalplane: 1 polygon for the whole footprint
        - ccd: 1 polygon per CCD (16 per footprint then)
        - quadrant: 1 polygon per quadrant (64 per per footprint then)
        ccd and quadrant level will account for gaps between the CCDs.
        
    Returns
    -------
    DataFrame
    """
    return FieldProjection.fieldid_to_radec(fieldid=fieldid, level=level, **kwargs)

def get_camera_location(ra, dec, level="quadrant", fielddf=None, keys=None):
    """ get expected location of targets withing the camera (focalplane, ccd, or quadrant) level.
    
    Parameters
    ----------
    ra, dec: float, array
        coordinates of interest (could be array of)

    level: str
        level of description of the camera.
        - focalplane: 1 polygon for the whole footprint
        - ccd: 1 polygon per CCD (16 per footprint then)
        - quadrant: 1 polygon per quadrant (64 per per footprint then)
        ccd and quadrant level will account for gaps between the CCDs.
        
    fielddf: pandas.DataFrame
        if available the results of 
        fielddf = ztffields.radec_to_fieldid(radec, level=level)
        Otherwise, it is just recomputed.

    keys: list, None
        name of the camera location.
        if None:
        - level='quadrant': [x,y] # pixels from the lower right corner of the quadrant
        - level='ccd': [i, j] # pixels from the lower right corner of the whole CCD
        - level='focalplane': [u, v] # in camera fraction

    Returns
    -------
    pandas.DataFrame
       - level quadrant: keys
    """
    radec = np.asarray([ra, dec]).T
    if fielddf is None:
        fielddf = ztffields.radec_to_fieldid(radec, level=level)
        
    
    coord_info = pandas.merge(fielddf, 
                                  pandas.DataFrame({"ra":ra, "dec":dec, "targetid":np.arange(len(dec))}), 
                                  left_index=True, right_index=True
                                 )
    # corner of fieldid
    lowleft_corner = FIELD_VERTICES.get_corner(level, which="lower left")
    upright_corner = FIELD_VERTICES.get_corner(level, which="upper right")
    
    fieldid_keys = lowleft_corner.index.names
    coord_info = coord_info.join(lowleft_corner, on=fieldid_keys, rsuffix="_lowleft"
                                 ).join(upright_corner, on=fieldid_keys, rsuffix="_upright")
    
    offset_percent = np.asarray([(coord_info["ra"]-coord_info["ra_lowleft"]
                                 )/(coord_info["ra_upright"]-coord_info["ra_lowleft"]),
                                 (coord_info["dec"]-coord_info["dec_lowleft"]
                                 )/(coord_info["dec_upright"]-coord_info["dec_lowleft"])])
    if level == "quadrant":
        if keys is None:
                keys = ["x", "y"]
        offset_percent *= np.asarray([RA_SIZE, DEC_SIZE])[:, None]    
    elif level == "ccd":
        if keys is None:
                keys = ["i", "j"]
        offset_percent *= np.asarray([RA_SIZE, DEC_SIZE])[:, None]*2
    elif level == "focalplane":
        if keys is None:
            keys = ["u", "v"]
        
    coord_info[keys] = offset_percent.T
    return coord_info.set_index(["targetid"]+fieldid_keys)[["ra","dec"]+keys]


# ============== #
#                #
# Generic Tools  #
#                #
# ============== #
def project_to_radec(verts_or_polygon, ra, dec):
    """ project a geometry (or its vertices) to given ra, dec coordinates
    
    Parameters
    ----------
    verts_or_polygon: shapely.Polygon or 2d-array
        geometry or vertices representing the camera footprint in the sky
        if vertices, the format is: x, y = vertices

    ra: float or array
        poiting(s) R.A.

    dec: float or array
        poiting(s) declination

    Returns
    -------
    list
        if input are vertices
        - list of new verticies
        if input are geometry
        - list of new geometries

    """
    from .utils import rot_xz_sph
    
    if type(verts_or_polygon) == geometry.Polygon: # polygon
        as_polygon = True
        fra, fdec = np.asarray(verts_or_polygon.exterior.xy)
    else:
        as_polygon = False
        fra, fdec = np.asarray(verts_or_polygon)
    
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)
    ra_, dec_ = np.squeeze(rot_xz_sph((fra/np.cos(fdec*np.pi/180))[:,None], 
                                            fdec[:,None], 
                                            dec)
                          )
    ra_ += ra
    pointings = np.asarray([ra_, dec_]).T
    if as_polygon:
        return [geometry.Polygon(p) for p in pointings]
    
    return pointings



def spatialjoin_radec_to_fields(radec, fields,
                                how="inner", predicate="intersects",
                                index_radec="index_radec",
                                allow_dask=True, **kwargs):
    """ join the radecs with the fields

    Parameters
    ----------
    radec: DataFrame or 2d-array 
        coordinates of the points. 
        - DataFrame: must have the "ra" and "dec" columns. 
            This will use the DataFrame's index are data index.
        - 2d array (shape N,2): returned index will be 'range(len(ra))'
    
    fields : [geopandas.geoserie, geopandas.geodataframe or  dict]
        fields contains the fieldid and fields shapes. Several forms are accepted:
        - dict: {fieldid: 2d-array, fieldid: 2d-array ...}
            here, the 2d-array are the field's vertices.

        - geoserie: geopandas.GeoSeries with index as fieldid and geometry as field's vertices.
            
        - geodataframe: geopandas.GeoDataFrame with the 'fieldid' column and geometry as field's vertices.

    Returns
    -------
    GeoDataFrame 
        (geometry.sjoin result)
    """
    # -------- #
    #  Coords  #
    # -------- #
    if type(radec) in [np.ndarray, list, tuple]:
        inshape = np.shape(radec)
        if inshape[-1] != 2:
            raise ValueError(f"shape of radec must be (N, 2), {inshape} given.")
        
        radec = pandas.DataFrame(np.atleast_2d(radec), columns=["ra","dec"])

    # Points to be considered
    geoarray = geopandas.points_from_xy(*radec[["ra","dec"]].values.T)
    geopoints = geopandas.GeoDataFrame({index_radec:radec.index}, geometry=geoarray)
    
    # -------- #
    # Fields   #
    # -------- #
    # goes from dict to geoseries (more natural) 
    fields = parse_fields(fields)

    # -------- #
    # Joining  #
    # -------- #
    # This goes linearly as size of fields
    if len(fields)>30_000 and allow_dask:
        try:
            import dask_geopandas
        except:
            pass # no more warnings.
            #warnings.warn("you don't have dask-geopandas, you have more than 20_000 fields, this would be faster with dask-greopandas")
        else:
            if type(fields.index) is not pandas.MultiIndex: # not supported
                fields = dask_geopandas.from_geopandas(fields, npartitions=10)
                geopoints = dask_geopandas.from_geopandas(geopoints, npartitions=10)
            else:
                warnings.warn("cannot use dask_geopandas with MultiIndex fields dataframe")
            
    sjoined = geopoints.sjoin(fields,  how="inner", predicate="intersects", **kwargs)
    if "dask" in str( type(sjoined) ):
        sjoined = sjoined.compute()

    # multi-index
    if type(fields.index) == pandas.MultiIndex:
        sjoined = sjoined.rename({f"index_right{i}":name for i,name in enumerate(fields.index.names)}, axis=1)
    else:
        sjoined = sjoined.rename({f"index_right": fields.index.name}, axis=1)

    return sjoined
    
def parse_fields(fields):
    """ read various formats for fields and returns it as a geodataframe

    Parameters
    ----------
    fields : [geopandas.geoserie, geopandas.geodataframe or  dict]
        fields contains the fieldid and fields shapes. Several forms are accepted:
        - dict: {fieldid: 2d-array or regions, fieldid: 2d-array or regions ...}
            here, the 2d-array are the field's vertices or a astropy/ds9 regions

        - geoserie: geopandas.GeoSeries with index as fieldid and geometry as field's vertices.
            
        - geodataframe: geopandas.GeoDataFrame with the 'fieldid' column and geometry as field's vertices.

    Returns
    -------
    GeoDataFrame (geometry.sjoin result)


    Examples
    --------
    provide a dict of ds9 regions
    >>> fields = {450:"box(50,30, 3,4,0)", 541:"ellipse(190,-10,1.5,1,50)"}
    >>> geodf = parse_fields(fields)

    """
    if type(fields) is dict:
        values = fields.values()
        indexes = fields.keys()
        # dict of array goes to shapely.Geometry as expected by geopandas
        test_kind = type( values.__iter__().__next__() ) # check the first
        if test_kind in [np.ndarray, list, tuple]:
            values = [geometry.Polygon(v) for v in values]
            
        if test_kind is str or "regions.shapes" in str(test_kind):
            values = [regions_to_shapely(v) for v in values]
            
        fields = geopandas.GeoSeries(values,  index = indexes)
            
    if type(fields) is geopandas.geoseries.GeoSeries:
        fields = geopandas.GeoDataFrame({"fieldid":fields.index},
                                        geometry=fields.values)
    elif type(fields) is not geopandas.geodataframe.GeoDataFrame:
        raise ValueError("cannot parse the format of the input 'fields' variable. Should be dict, GeoSeries or GeoPandas")

    return fields


def regions_to_shapely(region):
    """ 
    Parameters
    ----------
    region: str or Regions (see astropy-regions.readthedocs.io)
        if str, it is assumed to be the dr9 ircs format 
        e.g. region = box(40.00000000,50.00000000,5.00000000,4.00000000,0.00000000)
        if Regions, region will be converted into the str format
        using ``region = region.serialize("ds9").strip().split("\n")[-1]``
        The following format have been implemented:
        - box
        - circle
        - ellipse
        - polygon
        
    Returns
    -------
    Shapely's Geometry
        the geometry will depend on the input regions.
        
    Raises
    ------
    NotImplementedError
        if the format is not recognised.
        
    Examples
    --------
    >>> shapely_ellipse = regions_to_shapely('ellipse(54,43.4, 4, 2,-10)')
    >>> shapely_rotated_rectangle = regions_to_shapely('box(-30,0.4, 4, 2,80)')
    """
    import shapely
    from shapely import geometry
    
    if "regions.shapes" in str(type(region)):
        # Regions format -> dr9 icrs format
        region = region.serialize("ds9").strip().split("\n")[-1]

    tregion = type(region)
    if tregion is not str:
        raise ValueError(f"cannot parse the input region format ; {tregion} given")
        
    # it works, let's parse it.
    which, params = region.replace(")","").split("(")
    params = np.asarray(params.split(","), dtype="float")
    
    # Box, 
    if which == "box": # rectangle
        centerx, centery, width, height, angle = params
        minx, miny, maxx, maxy = centerx-width, centery-height, centerx+width, centery+height
        geom = geometry.box(minx, miny, maxx, maxy, ccw=True)
        if angle != 0:
            geom = shapely.affinity.rotate(geom, angle)
            
    # Cercle        
    elif which == "circle":
        centerx, centery, radius = params
        geom = geometry.Point(centerx, centery).buffer(radius)
        
    # Ellipse
    elif which == "ellipse":
        centerx, centery, a, b, theta = params
        # unity circle
        geom = geometry.Point(centerx, centery).buffer(1)
        geom = shapely.affinity.scale(geom, a,b)
        if theta != 0:
            geom = shapely.affinity.rotate(geom, theta)
        
    # Ellipse        
    elif which == "polygon":
        params = (params + 180) %360 - 180
        coords = params.reshape(int(len(params)/2),2)
        geom = geometry.Polygon(coords)
        
    else:
        raise NotImplementedError(f"the {which} form not implemented. box, circle, ellpse and polygon are.")
    
    # shapely's geometry
    return geom


    
class FieldProjection( object ):

    def __init__(self, level=None, radec=None):
        """ """
        self._fields = Fields(load_level=level)
        if radec is not None:
            self.load_projection(radec, level=level)

    # ============= #
    #  ClassMethod  #
    # ============= #
            
    @classmethod
    def radec_to_fieldid(cls, radec, level="focalplane", explode=True):
        """ """
        this = cls(radec=radec, level=level)
        return this.get_target_fields(explode=explode)
        
    @staticmethod
    def fieldid_to_radec(fieldid=None, level="focalplane", as_geoms=True):
        """ get centroid coordinated for given fieldid (and potentially CCD or quadrant sub-structures)
    
        Parameters
        ----------
        fieldid: int, list, array
            If None, all is returned.
            expected shape of fieldid depends on requested level:
            - focalplane: 1d array (N,) will return N entries e.g.
                (fieldid_1, fieldid_2, ...)
            - ccd, quadrant: 1 (N,) or 2d array (N, 2) will return N entries 
                if (N,): expected: (fieldid_1, 
                                    fieldid_2, ...)
                else: expected: ( [fieldid_1,ccdid_1], [fieldid_1,ccdid_2], 
                                  [fieldid_2,ccdid_3], ...)
                                  
            # Caution example:
            fieldid_to_radec([427,14], level="ccd") 
            => returns 32 entries (all 16 CCDs for fields 427 and 14
            fieldid_to_radec([ [427,14] ], level="ccd") 
            => returns 1 entry, the CCD 14 of field 427
            
            # Alternative
            for more option, simply use fieldid = None and parse the 
            pandas multi-index dataframe after. 
            For instance, to get all 13 of all fields:
            fieldid_to_radec(level="ccd").xs(13, level=1).
            
        level: str
            level of description of the camera.
            - focalplane: 1 polygon for the whole footprint
            - ccd: 1 polygon per CCD (16 per footprint then)
            - quadrant: 1 polygon per quadrant (64 per per footprint then)
            ccd and quadrant level will account for gaps between the CCDs.
            
        Returns
        -------
        DataFrame
        """
        which = FIELD_VERTICES.get_centroid(level, as_dataframe=True)
    
        if fieldid is None:
            return which
    
        fieldid = np.atleast_1d(fieldid)
        if fieldid.ndim == 1: # list of fields
            output = which.loc[fieldid] 
        else:
            if level == "focalplane":
                raise ValueError(f"{fieldid=} has {fieldid.ndim} dimensions, 1 expected for level='focalplane'.")
            elif fieldid.ndim == 2: # list of (fields, ccdid_rcid)
               output = which.loc[which.index.isin(fieldid.tolist())]
            else:
                raise ValueError(f"{fieldid=} has{fieldid.ndim} dimensions, 1 or 2 expected.")
            
        return output
        
    # ============= #
    #  Methods      #
    # ============= #        
    # - Top Level
    def get_target_fields(self, explode=False):
        """ """
        field_ = self.projection.groupby("index_radec")["fieldid"].apply(list)
        if explode:
            field_ = field_.explode()
            if self.level=="ccd":
                key = "ccdid"
            elif self.level=="quadrant":
                key = "rcid"
            else:
                key = None
            
            if key is not None:
                field_ = field_.str.split("_", expand=True).rename({0:"fieldid",1:key}, axis=1)

            field_ = field_.astype(int)
            
        return field_ 
    
    def get_field_targets(self, explode=False):
        """ """
        targets = self.projection.groupby("fieldid")["index_radec"].apply(list)
        if explode:
            targets = targets.reset_index().explode("index_radec")
            if self.level=="ccd":
                key = "ccdid"
                
            elif self.level=="quadrant":
                key = "rcid"
                
            else:
                key = None
            if key is not None:
                targets = targets.join(targets.pop("fieldid").str.split("_", expand=True).rename(
                                        {0:"fieldid",1:key}, axis=1)).set_index(["fieldid",key])
                targets = targets.astype(int)
        return targets
    
    # - mid Level  
    def load_projection(self, radec, level):
        """ """
        if type(level) != str:
            raise ValueError(f"level must be a str ; {type(level)} given")
        
        self._radec = radec
        self._projection = self._project_radec(radec, level=level)
        self._level = level
        
    def change_level(self, level):
        """ """
        if not hasattr(self, "_radec"):
            raise AttributeError("unknown radec, run load_projection first.")
            
        if level == self.level: # must have one otherwise no radec
            warnings.warn("no level to change")
            
        self.load_projection(self._radec, level=level)
        
    def get_geoseries(self, level):
        """ """
        geodf = self.fields.get_level_geodf(level)
        if type(geodf.index) == pandas.MultiIndex:
            for i in range(len(geodf.index.levshape)):
                level_ = geodf.index.get_level_values(i).astype(str)
                if i == 0:
                    index = level_
                else:
                    index +="_"+level_
                    
        else: # already single index
            index = geodf.index
        
        return geopandas.GeoSeries(geodf["geometry"].values, index=index)
        
        
    #
    # Core projection tool
    #
    def _project_radec(self, radec, level="ccd"):
        """ """
        geoserie = self.get_geoseries(level)
        return spatialjoin_radec_to_fields(radec, geoserie)
                
    # ============= #
    #  Properties   #
    # ============= #
    @property
    def fields(self):
        """ """
        return self._fields
    
    @property
    def projection(self):
        """ """
        if not hasattr(self,"_projection"):
            raise AttributeError("No projection loaded. See load_projection()")
        return self._projection
    
    @property
    def level(self):
        """ """
        if not hasattr(self,"_level"):
            raise AttributeError("No projection loaded. See load_projection()")        
            
        return self._level
    
    @property
    def target_radec(self):
        """ """
        if not hasattr(self,"_radec"):
            raise AttributeError("No projection loaded. See load_projection()")
            
        return self._radec
