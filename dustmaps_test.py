"""



"""

from __future__ import print_function
import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, QTable
import matplotlib.pyplot as plt
import scipy.stats as scistats
#import seaborn as sns
import astropy

#import passband_model_convolution as pmc
#import gaia_extinction
#import wdatmos

from dustmaps.bayestar import BayestarWebQuery as bwq

#bayestar = bwq(version= 'bayestar2017')
bayestar = bwq(version= 'bayestar2015')

ra = "14:31:44.6177"
dec= "-47:15:27.514"
#coords= coord.SkyCoord(ra, dec, unit= (u.hourangle, u.deg), distance = 1.5*u.kpc, frame = 'icrs')
coords= coord.SkyCoord(ra, dec, unit= (u.hourangle, u.deg), frame = 'icrs')

reddening, flags = bayestar(coords, mode= 'median', return_flags= True)
print(coords)
print(reddening)
print(flags['min_reliable_distmod'])
print(flags['max_reliable_distmod'])
