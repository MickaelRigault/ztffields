
import pandas
import numpy as np

from .fields import Fields

__all__ = ["skyplot_fields", "show_ztf_footprint"]

def values_to_color(value, cmap, vmin=None, vmax=None, alpha=None):
    """ """
    value = np.asarray(value)
    if vmin is None:
        vmin = np.nanmin(value)
    if vmax is None:
        vmax = np.nanmax(value)
        
    vrange = (value-vmin)/(vmax-vmin)
    return [vmin,vmax], vrange, cmap(vrange, alpha=alpha) 

def colorbar(ax, cmap, vmin=0, vmax=1, label="",
             fontsize="x-large", **kwargs):
    """ Set a colorbar in the given axis

    Parameters
    -----------
    ax: [mpl's Axes]
        Axis in which the colorbar will be drawn

    cmap: [mpl's colormap]
        A matplotlib colormap

    vmin, vmax: [float,float] -optional-
        Extend of the colormap, values of the upper and lower colors

    label, fontsize: [string, string/float] -optional-
        Label of the colorbar and its associated size
     
    **kwargs goes to matplotlib.colobar.ColorbarBase

    Return
    ------
    colorbar
    """
    import matplotlib
    
    if "orientation" not in kwargs.keys():
        bbox = ax.get_position()
        orientiation = "vertical" if bbox.xmax - bbox.xmin < bbox.ymax - bbox.ymin \
          else "horizontal"
        kwargs["orientation"] = orientiation

    norm    = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    c_bar   = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
                              norm=norm,**kwargs)
    
    c_bar.set_label(label,fontsize=fontsize)
    if "ticks" in kwargs.keys() and "ticklabels" not in kwargs.keys():
        c_bar.ax.set_xticklabels([r"%s"%v for v in kwargs["ticks"]])
        
    return c_bar



# -------------- #
#                #
#  FootPring     #
#                #
# -------------- #
def show_ztf_footprint(ax=None, ccdid=None, fieldid=None,
                           ccd_color="k", focalplane_color="0.7"):
    """ displays the ZTF footprint

    Parameters
    ----------
    ccdid: int, list, None
        if of the specific CCD(s) you want to display.
        None means all.

    fieldid: int
        provide a fieldid. The footprint will be 
        that of this field. If None, the camera 
        footprint will be at ra=0, dec=0.

    Returns
    -------
    figure
        
    """
    from .utils import ccdid_qid_to_rcid
    from matplotlib.patches import Polygon

    # consider as iterable for what comes next
    if ccdid is None:
        ccdids = np.arange(1,17)
    else:
        ccdids = np.atleast_1d(ccdid)

    if fieldid is None:
        ra, dec = 0,0
    else:
        ra, dec = Fields.get_field_centroids(fieldid).values
        
    # Figure
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure

    
    # Focal Plane
    fp = Fields.get_focalplane_contours(ra=ra, dec=dec)
    polyfp = Polygon(fp, facecolor="None", edgecolor=focalplane_color,
                         zorder=6, lw=1)
    ax.add_patch(polyfp)

    # Loops over the quadrant inside the CCDS
    # to make sure that they match as they should.
    # (it enabled a orientation debug, see get_layout())
    
    # CCD
    for ccdid in ccdids:
        ccd = Fields.get_ccd_contours(ccdid, ra=ra, dec=dec)
        polyccd = Polygon(ccd, facecolor="None", edgecolor="k", zorder=5, lw=2)
        ax.add_patch(polyccd)
        #centroid_ = 
        ax.text(*np.mean(ccd, axis=0), f"{ccdid}", va="center", ha="center", 
                color=ccd_color,zorder=5, weight="bold")
        
        # it's quadrants
        rcids = ccdid_qid_to_rcid(ccdid, np.arange(1,5))
        for i, rcid in enumerate(rcids):
            qad = Fields.get_quadrant_contours(rcid, ra=ra, dec=dec)
            polyqad = Polygon(qad, facecolor=f"C{i}", edgecolor="None", zorder=3, alpha=0.6)
            ax.add_patch(polyqad)
            centroid_ = np.mean(qad, axis=0)
            ax.text(*centroid_, f"{rcid}", va="center", ha="center", 
                    color="w", zorder=4, fontsize="x-small")


    ax.set_xlim(-5+ra,5+ra)
    ax.set_ylim(-5+dec,5+dec)
    return fig

# -------------- #
#                #
#  FieldFigure   #
#                #
# -------------- #
def skyplot_fields(fieldid, figsize=(7,4), level="focalplane",
                       system="icrs", projection=None, **kwargs):
    """ 
    
    Parameters
    ----------
    fieldid: list or pandas.Series
        index of the fields to show (see self.geodf).
        Two formats are accepted:
        - list: fieldid index, see facecolors for patches colors.
        - pandas.Series: fieldid will come from the index, the facecolor
        will be derived from the series' value.
        
    figsize : (float, float)
        Width, height in inches.

    level: str
        level of description of the camera.
        - focalplane: 1 polygon for the whole footprint
        - ccd: 1 polygon per CCD (16 per footprint then)
        - quadrant: 1 polygon per quadrant (64 per per footprint then)
        ccd and quadrant level will account for gaps between the CCDs.
    
    system: str
        coordinates system: 
        - icrs (ra, dec)
        - galactic (l, b)

    projection : matplotlib's string or cartopy.
        - if matplotlib:
        {None, 'aitoff', 'hammer', 'lambert', 'mollweide', 'polar', 'rectilinear', str}
        The projection type of the `~.axes.Axes`. *str* is the name of
        a custom projection, see `~matplotlib.projections`. The default
        None results in a 'rectilinear' projection.
        - if cartopy:
        any of the cartopy project from cartopy.crs.
            
    **kwargs goes to FieldFigure.skyplot_fields | e.g. add_mw, colorbar, edgecolor etc.

    Returns
    -------
    `matplotlib.figure`
    
    """
    return FieldFigure.skyplot_fields(fieldid, figsize=figsize, level=level,
                                      system=system, projection=projection,
                                      **kwargs)


class FieldFigure( object ):

    def __init__(self, geodf, figsize=(7,4), system="icrs", origin=180, **kwargs):
        """ 
        
        Parameters
        ----------

        geodf: geopandas.DataFrame
        

        figsize : (float, float)
            Width, height in inches.

        level: str
            level of description of the camera.
            - focalplane: 1 polygon for the whole footprint
            - ccd: 1 polygon per CCD (16 per footprint then)
            - quadrant: 1 polygon per quadrant (64 per per footprint then)
            ccd and quadrant level will account for gaps between the CCDs.

        system

        Returns
        -------
        None
        """
        self._geodf = geodf
        self.set_system(system, origin=origin)
        
        import matplotlib.pyplot as plt
        self._figure = plt.figure(figsize=(7,4))
        self._plotting = {}


    @classmethod
    def from_level(cls, level, figsize=(7,4), system="icrs", origin=180, **kwargs):
        """ """
        geodf = Fields.get_field_geometry(level=level, **kwargs)
        return cls(geodf=geodf, figsize=figsize, system=system, origin=origin)


    @staticmethod
    def _parse_system(system, default_frame="astro"):
        """ clean the input system for self consistancy 

        Parameters
        ----------
        system: str
            coordinates system: 
            - icrs (ra, dec) | shortcut for astro icrs
            - astro icrs (ra, dec) east is right
            - geo icrs (ra, dec) west is right
            - galactic (l, b)
            
            Note: for simplicity radec<->ircs or lb<->galactic
            such that 'geo radec' is equivalent to 'geo icrs'
            
        Returns
        -------
        system: str
        """
        system = system.replace("radec", "icrs").replace("lb", "galactic")
        *frame, sys = system.split(" ")
        if len(frame) == 0:
            frame = default_frame
        elif len(frame) == 1 and frame[0] in ["geo", "astro"]:
            frame = frame[0]
        else:
            raise ValueError(f"cannot parse input frame, astro or geo expected, got: {frame}")

        if sys == "galactic":
            return sys
        
        return f"{frame} {sys}"
        
        
    def set_system(self, system="icrs", origin=180):
        """ set the plotting coordinate system.

        Parameters
        ----------
        system: str
            coordinates system: 
            - icrs (ra, dec) | shortcut for astro icrs
            - astro icrs (ra, dec) east is right
            - geo icrs (ra, dec) west is right
            - galactic (l, b)

        origin: float
            where the center is (left is origin from 0)
            180 means that coordinates goes from -180->180
            
        Returns
        -------
        None
        """
        system = self._parse_system(system)
        
        xy_icrs = np.stack(self._geodf["geometry"].apply(lambda x: ((np.asarray(x.exterior.xy)).T) ).values)
        if "astro" in system:
            # invert ra direction
            xy_icrs[:,:,0] *=-1

        if "icrs" in system:
            xy = xy_icrs
            
        else:
            from astropy import coordinates, units            
            shape_ = xy_icrs.shape
            new_shape = np.prod(shape_[:-1]), shape_[-1]
            icrs = coordinates.ICRS( *xy_icrs.reshape(new_shape).T*units.degree )
            if system == "galactic":
                # Create a ICRS object and convert it to galactic coords
                gal = icrs.transform_to(coordinates.Galactic())
                # avoids edge effect.
                xy = np.asarray([gal.l.value, gal.b.value]).T.reshape(shape_)
                
            else:
                raise NotImplementedError(f"'{system}' has not been implemented")

        # identify and correct edge effects            
        self._geodf["xy"] = list(xy)
        self._system = system
        self._origin = origin
        
    # =============== #
    #   Class method  #
    # =============== #    
    @classmethod
    def skyplot_fields(cls, fieldid, figsize=(7,4), level="focalplane",
                           system="icrs", origin=180, projection=None,
                           add_mw=True, mwcolor = "0.7", mwb=5,
                           colorbar="hist", facealpha=0.9, edgecolor="None",
                           bins=None, **kwargs):
        """ plot fields on a 2d sky projection.

        This uses self.plot_sky()

        Parameters
        ----------
        fieldid: list or pandas.Series
            index of the fields to show (see self.geodf).
            Two formats are accepted:
            - list: fieldid index, see facecolors for patches colors.
            - pandas.Series: fieldid will come from the index, the facecolor
            will be derived from the series' value.
        
        figsize : (float, float)
            Width, height in inches.

        level: str
            level of description of the camera.
            - focalplane: 1 polygon for the whole footprint
            - ccd: 1 polygon per CCD (16 per footprint then)
            - quadrant: 1 polygon per quadrant (64 per per footprint then)
            ccd and quadrant level will account for gaps between the CCDs.
    
        system: str
            coordinates system: 
            - icrs (ra, dec)
            - galactic (l, b)

        projection : matplotlib's string or cartopy.
            - if matplotlib:
            {None, 'aitoff', 'hammer', 'lambert', 'mollweide', 'polar', 'rectilinear', str}
            The projection type of the `~.axes.Axes`. *str* is the name of
            a custom projection, see `~matplotlib.projections`. The default
            None results in a 'rectilinear' projection.
            - if cartopy:
            any of the cartopy project from cartopy.crs.

        add_mw: bool
            should the Milky Way location be indicated

        mwcolor: matplotlib's color
            = ignored if add_mw=False =
            color of the Milky Way

        mwb: float
            = ignored if add_mw=False =
            width of the Milky Way (in degree from b=0)

        colorbar: str or bool
            should this add a colorbar and which ?
            - colorbar='hist': histcolorbar
            - else, if not None: colorbar

        edgecolor: matplotlib's color
            edgcolor of the patches. "None" means no color

        facealpha: float
            alpha value for the patches.

        bins: int
            sequencing of the colormap. 'lut' in matplotlib's get_cmap()

        **kwargs goes to plot_sky | e.g. cmap, vmin, vmax ...

        Returns
        -------
        figure

        """
        # Create figure and the plotting system        
        figfield = cls.from_level(level=level, figsize=figsize, system=system, origin=origin)
        
        # Create the axes
        ax = figfield.add_axes(projection=projection)
        
        # Draw the fields as matplolib's polygon
        _ = figfield.plot_sky(ax=ax, fieldid=fieldid,
                              facealpha=facealpha, edgecolor=edgecolor,
                              colorbar=colorbar, bins=bins, **kwargs)
        
        # Add the Milky way if wanted.
        if add_mw:
            _ = figfield.add_milkyway(ax=ax, color=mwcolor, b=mwb)

        return figfield.figure
    
    # =============== #
    #   add other     #
    # =============== #
    def add_axes(self, rect=[0.15,0.22,0.75,0.75],
                     projection=None, 
                     reset=False, **kwargs):
        """ add an axes
        
        Parameters
        ----------
        rect : sequence of float
            The dimensions [left, bottom, width, height] of the new Axes. All
            quantities are in fractions of figure width and height.

        projection : matplotlib's string or cartopy.
            - if None:
               if cartopy is installed this will use mollweide from cartopy
               else from matplotlib.
            - if matplotlib:
            {None, 'aitoff', 'hammer', 'lambert', 'mollweide', 'polar', 'rectilinear', str}
            The projection type of the `~.axes.Axes`. *str* is the name of
            a custom projection, see `~matplotlib.projections`. The default
            None results in a 'rectilinear' projection.
            - if cartopy:
            any of the cartopy project from cartopy.crs.

        reset: bool
            Should this reset the current plotting ?
            
        **kwargs goes to matplotlib's figure.add_axes()
        
        Returns
        -------
        `~.axes.Axes`, or a subclass of `~.axes.Axes`
            The returned axes class depends on the projection used. It is
            `~.axes.Axes` if rectilinear projection is used and
            `.projections.polar.PolarAxes` if polar projection is used.
            
        """
        if reset:
            self.reset()
            
        if not hasattr(self,"_figure") or self.figure is None:
            import matplotlib.pyplot as plt
            self._figure = plt.figure(figsize=(7,4))

        if projection is None:
            try:
                import cartopy.crs as ccrs
                projection = ccrs.Mollweide()
            except ImportError:
                projection = "mollweide"
            
        # check if input projection is from cartopy
        self._is_cartopy = "cartopy" in str( type(projection) )
        
        ax = self.figure.add_axes(rect, projection=projection, **kwargs)

        if self._is_cartopy:
            self._projection = projection
            ax.set_global()
#            self.set_system(self._system)
            
        return ax
    
    # =============== #
    #     Main        #
    # =============== #        
    def plot_sky(self, fieldid=None, ax=None, 
                 facecolors=None, system=None,
                 cmap="cividis", bins=None, vmin=None, vmax=None,
                 facealpha=None, colorbar=None, label=None,
                 **kwargs):
        """ Plot fields on a 2d sky projection
        
        Parameters
        ----------
        fieldid: list or pandas.Series
            index of the fields to show (see self.geodf).
            Two formats are accepted:
            - list: fieldid index, see facecolors for patches colors.
            - pandas.Series: fieldid will come from the index, the facecolor
            will be derived from the series' value.
            
        ax: matplotlib.axes
            axes where figure will be shown. If None, this call
            self.add_axes()
            
        facecolors: array
            = ignored if fieldid is a pandas.Series = 
            value used to build the patches color. 
            len of facecolors must match len of fieldid.
        
        system: str
            coordinates system: 
            - icrs (ra, dec)
            - galactic (l, b)

        cmap: matplotlib.Colormap or str
            if str: a name of a colormap known to Matplotlib.
            
        bins: int
            sequencing of the colormap. 'lut' in matplotlib's get_cmap()
        
        vmin: float
            minimal value for the colormap
            
        vmax: float
            maximal value for the colormap
            
        facealpha: float
            alpha value for the patches.
            
        edgecolor: matplotlib's color
            edgcolor of the patches. "None" means no color

        colorbar: str or bool
            should this add a colorbar and which ?
            - colorbar='hist': histcolorbar
            - else, if not None: colorbar

        label: str
            label of the colorbar

        Returns
        -------
        figure
            equivalent to self.figure
            
        Example
        -------
        _ = figfield.plot_sky(fieldid=fieldid_s, 
                              facealpha=0.9, edgecolor="None",
                              colorbar="hist")

        """
        if type(fieldid) == pandas.Series:
            facecolors = fieldid.values
            fieldid = fieldid.index

        if system is not None:
            self.set_system(system)
            
        import matplotlib.pyplot as plt

        fieldverts = self.geodf[["xy"]].copy()
        if fieldid is not None:
            fieldverts = fieldverts.loc[fieldid]
            
        if ax is None:
            ax = self.add_axes(reset=True)

        self._plotting["ax"] = ax
            
        #
        # Per-Polygon options
        #
        # - facecolors
        prop = {} # generic properties applicable to all fields
        if facecolors is not None and len(facecolors) == len(fieldverts):
            cmap = plt.get_cmap(cmap, lut=bins)
            vrange, cdata, colors = values_to_color(facecolors, 
                                        cmap=cmap, 
                                        vmin=vmin, vmax=vmax,
                                        alpha=facealpha)
            self._plotting["bins"] = bins
            self._plotting["cmap"] = cmap
            self._plotting["vrange"] = vrange
            self._plotting["cdata"] = cdata
            fieldverts = fieldverts.join(pandas.DataFrame({"facecolor":list(colors)}, 
                                                            columns=['facecolor'], 
                                                            index=fieldverts.index))
            
        else:
            from matplotlib.colors import to_rgba
            self._plotting["bins"] = None
            prop["facecolor"] = to_rgba("C0",0.1)
            prop["edgecolor"] = to_rgba("k",0.7)
            
        #
        # Calling plotter
        #
        _ = self._show_fieldverts(fieldverts, ax=ax, **{**prop,**kwargs})
        if colorbar is not None:
            if type(colorbar) is str and "hist" in colorbar:
                cbar,_ = self.add_histcolorbar()
            else:
                cbar = self.add_colorbar()
            if label is not None:
                cbar.set_label(label, fontsize="large")
            
        return self.ax.figure
    
    # =============== #
    #   add other     #
    # =============== #    
    def add_colorbar(self, cax=None, ymin=0.15, height=0.04, **kwargs):
        """ add a colorbar
        
        This call colobar()
        
        Parameters
        ----------
        cax: matplotlib.axes
            axes where the colorbar should be displayed
            
        ymin: float
            = if cax is None =
            bottom coordinate of the axes to be created (cax)
            
        height: float
            = if cax is None =
            height coordinate of the axes to be created (cax)
            
        **kwargs goes to colorbar() | label, fontsize...
        
        Returns
        -------
        matplotlib.colorbar.ColorbarBase
            colorbar
        """
        if cax is None:
            bboxax = self.ax.get_position()
            fig = self.ax.figure
            cax = fig.add_axes([bboxax.xmin, ymin, bboxax.width, height])
            
        vmin, vmax = self._plotting["vrange"]
        return colorbar(cax, self._plotting["cmap"], 
                        vmin=vmin, vmax=vmax, **kwargs)
        
    def add_histcolorbar(self, cax=None, bins=None,
                            ymin=0.15, height=0.06, hratio=0.7, gapratio=0.05,
                            histalpha=None, **kwargs):
        """ add a colorbar with histogram on top.
        
        This call colobar()
        
        Parameters
        ----------
        cax: [matplotlib.axes,matplotlib.axes]
            axes where the colorbar AND the histogram should be displayed.
            `cax_, histax = cax`
            where cax_ is the actual colorbar and histax the axes for the
            histogram.
                                
        bins : int or str
            The number of equal-width bins in the given range.

        ymin: float
            = if cax is None =
            bottom coordinate of the axes to be created (cax)
                            
        height: float
            = if cax is None =
            height coordinate of the axes to be created (cax)        
            
        hratio: float
            = if cax is None =
            fraction of the height that goes to histogram axes
            
        gapratio: float
            = if cax is None =
            fraction of the height to skip (span) between the two axes.
            
        histalpha: float
            matplotlib's alpha on the histogram.
            
        **kwargs goes to colorbar() | label, fontsize...
        
        Returns
        -------
        colorbar, axes.bar's output.
        """
        import matplotlib.pyplot as plt
        if cax is not None:
            cax, histax = cax
        else:
            bboxax = self.ax.get_position()
            cax = self.ax.figure.add_axes([bboxax.xmin, ymin, bboxax.width, height*(1-(hratio+gapratio))])
            histax = self.ax.figure.add_axes([bboxax.xmin, ymin+height*(1-hratio), bboxax.width, 
                                         height*hratio])

        if "vrange" not in self._plotting:
            cax.set_visible(False)
            histax.set_visible(False)
            return None, None
    
        cax.set_visible(True) # just to make sure
        histax.set_visible(True) # just to make sure
        
        vmin, vmax = self._plotting["vrange"]
        data = self._plotting["cdata"]
        cdata = self.get_current_colordata()
        
        # build histogram
        if bins is None:
            bins = self._plotting.get("bins", "auto")
            if bins is None:
                bins = 'auto'

        self._build_histrogram(cdata, bins=bins)
        hbins = self._plotting["hist"]["bins"]
        sequence = hbins["size"]
        cmap = plt.get_cmap(self._plotting["cmap"].name, sequence)
        
        # Plotting colorbar
        cbar = colorbar(cax, cmap, vmin=vmin, vmax=vmax, **kwargs)
        
        hcolors = cmap((hbins["centroid"] - hbins["vmin"])/(hbins["vmax"] - hbins["vmin"]))
        _ = histax.bar(hbins["centroid"], self._plotting["hist"]["intensity"], 
                           width=hbins["width"], color=hcolors,
                           alpha=histalpha)
        histax.set_ylim(bottom=0)
        histax.set_xlim(vmin, vmax)
        histax.axis("off")
        return cbar, _
        
    def add_milkyway(self, b=5, nbins=400, l_start=0, l_stop=360, ax=None, zorder=2, **kwargs):
        """ add the milky way location on the main axes' figure (self.ax)
        
        Parameters
        ----------
        b: int
            width (in degree) of the band around b=0
            
        nbins: int
            number of bins one the b=0 line (resolution)
            
        l_start, l_stop: float
            l coordinate where milky way is defined.
            
        **kwargs goes to ax.plot (if b is None) or to ax.fill_between (otherwise)
        
        Returns
        -------
        output of ax.plot or ax.fill_between
        """
        from astropy import coordinates, units
        if ax is None:
            ax = self.ax

        
        prop_line = dict(ls="-", color="0.7", alpha=0.5)
        prop_band = dict(facecolor="0.7", alpha=1)
        kwargs["zorder"] = zorder # because explicit input

        if b is None or b==0:
                b = 1
        
        if "icrs" in self.system:
            long_, lat_ = np.linspace(l_start, l_stop,nbins)*units.deg, np.zeros(nbins)*units.deg
            gal = coordinates.Galactic( long_, lat_ ).transform_to(coordinates.ICRS())
            ra, dec = gal.ra.value, gal.dec.value
            if "astro" in self.system:
                ra *=-1
            # make sure nothign breaks

        ra, dec = self._parse_radec_for_plot(ra, dec, sort=True)        
        #              #
        # Cartopy      #
        #              #

        if self._is_cartopy:
                
            from shapely.geometry import LineString
            import cartopy.crs as ccrs
            
            if "icrs" in self.system:
                ggal = LineString( zip(ra, dec) ).buffer(b)
                
            elif self.system == "galactic":                
                ggal = LineString( zip(long_.value, lat_.value) ).buffer(b)
                
            else:
                raise NotImplementedError("Milky way plotting has not been implemented for the current system '{self.system}'")

            return ax.add_geometries(ggal,
                                     ccrs.PlateCarree(self._origin),
                                     **{**prop_band, **kwargs}
                                     )
            
        #              #
        # Matplotlib   #
        #              #
        # l, b
        if "galactic" in self.system:
            if b is None:
                _ = ax.axhline(0, **{**prop_line, **kwargs})
            else:
                _ = ax.axhspan(-b*np.pi/180, b*np.pi/180, **{**prop_band,**kwargs})
                
        # RA Dec
        elif "icrs" in self.system:
            # matplotlib works in radian
            _ = ax.fill_between(ra,
                                dec + b*np.pi/180,
                                dec - b*np.pi/180,
                                **{**prop_band, **kwargs})
                
        else:
            raise NotImplementedError("Milky way plotting has not been implemented for the current system '{self.system}'")
        
        return _

    def _parse_radec_for_plot(self, ra, dec, sort=False):
        """ """
        if self._is_cartopy:
            pass
        else:
            ra = ra+self._origin # copy
            ra = (ra + 180)%360 - 180 # -180->+180
            # copy
            ra = ra * np.pi/180
            dec = dec * np.pi/180

        if sort:
            argsort_ra = np.argsort(ra)
            ra = ra[argsort_ra]
            dec = dec[argsort_ra]

        return ra, dec
        
    
    # =============== #
    #   Internal      #
    # =============== #
    def get_current_colordata(self):
        """ """
        cdata = self._plotting.get("cdata")
        if cdata is None:
            return None
        vmin, vmax = self._plotting["vrange"]
        vscale = vmax-vmin
        return  cdata*vscale + vmin
        
    def reset(self):
        """ """
        self._plotting = {}
    
    def _show_fieldverts(self, fieldverts, ax, zorder=4, **kwargs):
        """ """
        kwargs["zorder"] = zorder # because explicit input
        
        from matplotlib.patches import Polygon
        if self._is_cartopy:
            import cartopy.crs as ccrs
            kwargs["transform"] = ccrs.PlateCarree(self._origin)

        if not self._is_cartopy:
            import warnings
            warnings.warn("you are not using cartopy projection, issues may appear on the edge")
            fieldverts = fieldverts.copy()
            xy = np.stack(fieldverts["xy"].values)
            xy[:,:, 0] += self._origin
            xy[:,:, 0] = (xy[:,:, 0] + 180)%360 - 180 # -180->+180
            xy *= np.pi/180
            fieldverts["xy"] = list(xy) # feed back in
            
        patches = [Polygon(**fieldvert_.to_dict(), **kwargs) for i, fieldvert_ in fieldverts.iterrows()]
        return [ax.add_patch(patch_) for patch_ in patches]
            

    def _build_histrogram(self, data, bins=None, **kwargs):
        """ """
        intensity, binegdes = np.histogram(data, range=self._plotting["vrange"], bins=bins,  **kwargs)
        bins = {"edge":binegdes,
                "centroid":np.mean([binegdes[1:],binegdes[:-1]], axis=0),
                "width":binegdes[1:]-binegdes[:-1],
                "size":len(binegdes)-1,
                "vmin":binegdes[0],
                "vmax":binegdes[-1]}
        
        self._plotting["hist"] = {"intensity": intensity,
                                  "bins": bins}
        
    # ============== #
    #   Properties   #
    # ============== #
    @property
    def geodf(self):
        """ """
        return self._geodf
    
    @property
    def figure(self):
        """ """
        return self._figure
    
    @property
    def ax(self):
        """ """
        return self._plotting.get("ax", None)
                              

    @property
    def system(self):
        """ coordinate system (icrs, galactic etc.) """
        return self._system
