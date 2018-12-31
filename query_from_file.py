"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-12-31

This should specifically be capable of taking in  a list of coordinates and such and then querying the gaia archive for the associated objects... hopefully. We'll see! This has only failed every other time I've attempted it.

"""



import numpy as np
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
 
