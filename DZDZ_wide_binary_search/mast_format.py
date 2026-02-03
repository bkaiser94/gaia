"""
Created by Ben Kaiser (UNC-Chapel Hill) 2026-01-29

Export the name and coordinates of the relevant objects as a single list that can be searched in MAST



"""

import numpy as np
#from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column, join
import time
import sys
start = time.time()

input_name='WDWD_widebinaries_eDR3_simbadadded1_simbadadded2_DZcands_Vincent2023bClass1_Vincent2023bClass2.csv'

input_table=Table.read(input_name)


name_list=[]
ra_list=[]
dec_list=[]
num_list=['1','2']
for num in num_list:
    ra_list.append(input_table['ra'+num])
    dec_list.append(input_table['dec'+num])
    name_list.append(input_table['simbad_name'+num])

ra_array=np.hstack(ra_list)
dec_array=np.hstack(dec_list)
name_array=np.hstack(name_list)
output_table=Table(names=['Target','ra','dec'],data=[name_array,ra_array,dec_array])
output_table.pprint

output_table.write('mast_formatted_DZcands.csv')
