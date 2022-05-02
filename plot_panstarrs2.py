"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-05-02.

Hopefully this will plot the panstarrs data for the purple objects that I retrieved from MAST.

Probably is going to have to be run in python3

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
#import astropy


#import passband_model_convolution as pmc
import gaia_extinction
#import wdatmos
import plotting_dicts as pod


