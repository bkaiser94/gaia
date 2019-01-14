"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-01-09

Should take one of the fuller Gaia tables (meaning like all of the columns for the small target numbers), and
output a SOAR-compatible target list (so space-delimited or tab)

"""


import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column

