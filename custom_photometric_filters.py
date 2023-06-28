"""
Created by Ben Kaiser (UNC-Chapel Hill) 2023-06-06 (D-Day).

This script is either going to create CSV files of concocted photometric systems or it is going to
construct them to be used in xp_plotting.py. These photometric systems will really work more
like spectroscopic indices because I will have them target series of absorption features.

There should be one set targeting the strong metals in low-Teff DZs, which would be Ca II H&K, 
Ca I resonance, and Na I D. There should be another set that targets the Balmer lines probably 
as a proof of concept that this kind of works too. Then there should be the third set that targets 
continuum around the other two sets. This set of 'filters' should then be implemented in 
xp_plotting.py as the custom photometric filters and combined with XP spectra.


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
import gaiaxpy as xpy
import pandas as pd


#import passband_model_convolution as pmc
#import gaia_extinction
#import wdatmos
import plotting_dicts as pod

#phot_system=xpy.PhotometricSystem.Gaia_DR3_Vega

#print(phot_system.__init__)

metal_band=[
    [3914.,3987.],
    [4179.,4339.],
    [5865.,5932.]
    ]

balmer_band=[
    [3827.,3857.],
    [3874.,3916.],
    [3952.,4000.],
    [4074.,4142.],
    [4293.,4382.],
    [4816.,4917.],
    [6511.,6610.]
    ]

continuum_band=[
    [3000.,3760.],
    [4500.,4600.],
    [5000.,5700.],
    [6000.,6450.],
    [6750.,10000.]
    ]

continuum_blue_band=[
    [3000.,3760.],
    [4500.,4600.],
    [5000.,5700.],
    [6000.,6450.]
    ]

continuum_red_band=[
    [6750.,10000.]
    ]


filter_dict={
    'metal':metal_band,
    'balmer':balmer_band,
    'continuum':continuum_band,
    'blue_continuum':continuum_blue_band,
    'red_continuum':continuum_red_band
    }

#print(metal_band)
for name in filter_dict:
    print(name,filter_dict[name])





















