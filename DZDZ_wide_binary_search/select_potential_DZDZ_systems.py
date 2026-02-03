"""
Created by Ben Kaiser (UNC-Chapel Hill) 2026-01-29.

This should identify all of the WD-WD binaries with at least one component that has a "Z" in its spectral type and then create a sublist of those to output. Then write that as an output.






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
from astropy.table import Table, Column, vstack, join, MaskedColumn
#import scipy.interpolate as scinterp
import time
#import pyvo
#from astroquery.vizier import Vizier

#from astroquery.gaia import Gaia
#from astroquery.xmatch import XMatch


sys.path.append('../')
sys.path.append('/Users/BenKaiser/Desktop/radial_velocity_calculations/')
import cal_params as cp

input_name='WDWD_widebinaries_eDR3_simbadadded1_simbadadded2.csv'
input_table=Table.read(input_name)

output_name=input_name.split('.')[0]+'_DZcands.csv'

num_list=['1','2']
spname='simbad_sp_type'

Z_list=[]
for num in num_list:
    safe=input_table[spname+num].filled('')
    mask=np.char.find(safe.astype(str),'Z')>=0
    filtered=input_table[mask]
    Z_list.append(filtered)
    
output_table=vstack(Z_list)

output_table.pprint()
print(output_table['simbad_sp_type1'],output_table['simbad_sp_type2'])

output_table.write(output_name)




