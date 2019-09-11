"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-09-06


Take a Gaia dataset (that includes all of the bells and whistles) and generate galactic kinematic quantities if they
aren't known and save them into the file in addition to the previously available values. If they are already known,
just go ahead and make the plots that show memberships and stuff.

This is probably not going to work correctly for awhile.

THIS SHOULD BE RUN IN PYTHON 3 BECAUSE ASTROPY NEEDS TO BE ON ITS GAME!!!

"""

from __future__ import print_function
import numpy as np
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, QTable
import matplotlib.pyplot as plt
import scipy.stats as scistats


import plotting_dicts as pod


print('\n\n YOU BETTER BE USING PYTHON3!\n\n\n')

list_color = '#1ca1f2'


#######3error distribution variables
mc_number = 10000
percent_off = 34 #1-sigma equivalent
#############

#target_input ='20190516B_retargeted_purple_search_gaia_scbd.csv'
target_input='20190829_alkaliWD_targeted_gaia_scbd.csv'

target_table = Table.read(target_input)




def get_galactic_coords(row):
    star_coord= coord.SkyCoord(row['ra']*u.deg, row['dec']*u.deg, pm_ra_cosdec= row['pmra']*u.mas/u.yr, pm_dec= row['pmdec']*u.mas/u.yr, frame='icrs')
    print('l', row['l'], 'b', row['b'])
    print(star_coord.galactic)
    print(star_coord.galactic.l.value - row['l'], star_coord.galactic.b.value - row['b'])
    
for row in target_table:
    get_galactic_coords(row)

