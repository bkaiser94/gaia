"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-06-21

This needs to be run in the python3.13 for me. ESOAsg requires python 3.7 at minimum to be installed.



"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
from astropy.table import Table, Column, vstack, join
#import scipy.interpolate as scinterp
import time
import pyvo
from astroquery.vizier import Vizier
from astroquery.eso import Eso

from astroquery.gaia import Gaia
from ESOAsg import archive_catalogues

#import requests

#url='http://archive.eso.org/tap_cat/async/phase'
#response = requests.get(url)
#response.raise_for_status()  # Raises an exception for 4xx and 5xx status codes
#data = response.json()  #

sys.path.append('../')


input_file='topcat_alleDR3_GentileFusillo__WDs_vagueVVVfootprint_Xmatch_VVVViracAstrometry_roughcorrPMRAroughcorrPMRDEsubset.fits'
input_directory='/Users/BenKaiser/Desktop/gaia/IR_variable_search/'

input_table=Table.read(input_directory+input_file)
input_table.pprint()

catalogues=archive_catalogues.catalogues_info(collections=['VVVX'])
catalogues['collection','title','version','description','table_name'].pprint()

for row in catalogues:
    print(row['table_name'])
print('Getting catalogue')
#vvvx_lc=archive_catalogues.get_catalogues(tables=['VVVX_VIRAC_V2_LC'])
vvvx_lc=archive_catalogues.get_catalogues(tables=['VVVX_VIRAC_V2_SOURCES'])
print('Got Catalogue')

vvvx_lc.pprint()

for colname in vvvx_lc.colnames:
    print(colname)
    
    
new_table=Table.join(input_table,vvvx_lc,keys_left='srcid',keys_right='sourceid')

new_table.pprint()
