"""
Created by Ben Kaiser (UNC-Chapel Hill) on 2024-10-16.

This script should take the ASCII tab-delimited list of objects from Simbad that I
output from the Kesseli et al. 2019 paper and format it in a way that I can open it in my Gaia querying. 


"""



import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column
import matplotlib.pyplot as plt

import time
import sys
start = time.time()




#input_file='Kesseli_2019_subdwarfs_Simbad_results.txt'
input_file='Kesseli_2019_subdwarfs_Simbad_results_otherheader.txt'

output_filename='Kesseli_2019_subdwarfs_J2016coords_simbad.csv'

input_table=Table.read(input_file, format='ascii.fixed_width', delimiter='\t')


input_table.pprint()

print(input_table.colnames)

colnames_all=[
    
    ]

good_coords=np.where(input_table['coord2 (ICRS,J2016/2000)']!='No Coord.')
sub_table=input_table[good_coords]

print(sub_table['coord1 (ICRS,J2000/2000)'])

coord_col=sub_table['coord2 (ICRS,J2016/2000)']

name_col=sub_table['identifier']
name_col.name='name'

spectype_col=sub_table['spec. type']
spectype_col.name='spec_type'

G_col=sub_table['Mag G']
G_col.name='phot_g_mean_mag'

ra_list=[]
dec_list=[]

for row in coord_col:
    print('row',row)
    parts=row.split(' ')
    print(parts)
    ra=':'.join(parts[:3])
    dec=':'.join(parts[3:])
    ra_list.append(ra)
    dec_list.append(dec)
    print('ra',ra)
    print('dec',dec)

ra_col=Column(ra_list, name='ra')
dec_col=Column(dec_list, name='dec')

#print(ra_col)
#print(dec_col)
#good_coords=np.where(coord_col!='No Coord.')

#bad_coords=np.where(coord_col=='No Coord.')

#print('bad_coords',bad_coords)

output_table=Table([name_col, ra_col, dec_col, G_col, spectype_col])
output_table.pprint()

print('Saving file:', output_filename)
output_table.write(output_filename, format='ascii.csv',overwrite=False)
print('File saved:', output_filename)
