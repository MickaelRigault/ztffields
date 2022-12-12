"""
projection
----------

Project targets (ra, dec) into the fields.

.. autosummary::

   radec_to_fieldid
   fieldid_to_radec
   spatialjoin_radec_to_fields
   FieldProjection

"""


import warnings
import pandas
import numpy as np
#
import geopandas
#
from .fields import Fields

__all__ = ["radec_to_fieldid", "fieldid_to_radec", "FieldProjection"]

# ============== #
#                #
#  shortcuts     #
#                #
# ============== #

def radec_to_fieldid(radecs, level="focalplane", **kwargs):
    """ """
    return FieldProjection.radec_to_fieldid(radecs, level=level, **kwargs)

def fieldid_to_radec(radecs, level="focalplane", **kwargs):
    """ """
    return FieldProjection.fieldid_to_radec(radecs, level=level, **kwargs)



# ============== #
#                #
# Generic Tools  #
#                #
# ============== #

def spatialjoin_radec_to_fields(radec, fields,
                                how="inner", predicate="intersects",
                                index_radec="index_radec", **kwargs):
    """ 

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
        if (inshape:=np.shape(radec))[-1] != 2:
            raise ValueError(f"shape of radec must be (N, 2), {inshape} given.")
        
        radec = pandas.DataFrame(radec, columns=["ra","dec"])

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
    return geopoints.sjoin(fields,  how="inner", predicate="intersects", **kwargs)


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
        
    @classmethod
    def fieldid_to_radec(cls, radec, level="focalplane", explode=True):
        """ """
        this = cls(radec=radec, level=level)
        return this.get_field_targets(explode=explode)
        
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
