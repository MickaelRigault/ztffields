
========================
ztffields documentation
========================

ztffields_ is a package made to handle Zwicky Transient Facility (ZTF)
fields system.

ZTF observes using pre-define sets of field idenfitied
with a unique `fieldid` ; it already exists more than 2000 fieldids.

The ZTF camera is made of 16 CCDs identified by a unique `ccdid`
ranging from 1 to 16 (incl.) and each CCD has four amplifier,
definining 4 quadrants defining, in total 64 rcid (16x4) ranging from
0 to 63 (incl.).


ztffields_ enables to interact with the ztf field system at the
focalplane (whole camera footprint ; 1 polygon), at the ccd (16 polygons) or are the
quadrant (64 polygons) levels.

Sharp start
============

**Get the fieldid and quadrant rcid that contains a given set of targets**

.. code-block:: python
		
	np.random.seed(1234)
	# generate 4000 RA, Dec coordinates.
	radecs = pandas.DataFrame({"ra":np.random.uniform(-30, 90, size=4000),
                           "dec":np.random.uniform(0, 150, size=4000)}
                         )
	df = ztffields.radec_to_fieldid(radecs, level="ccd") # ccd level
	df.loc[600] # see the fieldid and ccdid containing the target 600

Matching 4000 targets takes: 300ms (focalplane level), ~1s (ccd-level)
or ~4s (quadrant level)

.. toctree::
   :caption: Code documentation
   :titlesonly:

   ztffields
   ztffields.projection
   ztffields.fields
   ztffields.utils   
   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
