"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-01-09

Should take one of the fuller Gaia tables (meaning like all of the columns for the small target numbers), and
output a SOAR-compatible target list (so space-delimited or tab)

"""

from __future__ import print_function
import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column

input_file= '20190107_chris_merge_gaia.csv'
input_file='20190405_purple_search_gaia_unique.csv'
input_file='exc1_8_2_2_purple_search_gmaglimit_gaia_sc.csv'
#input_file='Lindegren_odd_survivors_gaia_sc.csv'
#input_file='sdssj1330_similar_subset_gaia_sc.csv'
#output_file= '20190107_chris_merge_targlist.txt'
#input_file='expanded_purple_search_gmaglimit_gaia_sc.csv'
comment_string='coolWDdMcand'

output_file= input_file.split('.')[0]+'_targlist.txt'
input_table= Table.read(input_file, format= 'ascii.csv')
delimiter='\t'

sort_by_ra= True

try:
    name_array= np.array(input_table['name'])
except KeyError as error:
    print(error)
    print('No names column, so just making all Gaia names')
    name_array=np.full(input_table['ra'].shape, '.')

coords= coord.SkyCoord(input_table['ra'], input_table['dec'], unit=(u.deg, u.deg))
string_coords= coords.to_string(style='hmsdms')
replace_chars= ['d','h','m']
ra_list=[]
dec_list=[]
epoch_list=[]
mag_list=[]
name_list=[]
mag_string='g'
full_mag_string= 'phot_'+mag_string+'_mean_mag'
for thing, name, mag  in zip(string_coords, name_array, input_table[full_mag_string]):
    print(name)
    print(thing)
    for char in replace_chars:
        thing= thing.replace(char, ':')
    thing=thing.replace('s', '')
    split_string=thing.split(' ')
    ra= split_string[0][:11] #limiting precision of the decimal
    dec= split_string[1][:12] #limiting precision of the decimal, need the additional index for sign
    small_ra= ra.replace(':', '')[:4]
    small_dec=dec.replace(':','')[:5] #need additional index for + or -
    if name=='.':
        name='GaiaJ'+small_ra+small_dec
        print('new name:', name)
    elif name=='none':
        name='GaiaJ'+small_ra+small_dec
        print('new name:', name)
    elif name == '0':
        name='GaiaJ'+small_ra+small_dec
        print('new name:', name)
    print(thing)
    name_list.append(name)
    ra_list.append(ra)
    dec_list.append(dec)
    epoch_list.append('2000.0')
    #mag_list.append(mag_string.upper()+'='+str(mag)[:4])
    mag_list.append(mag_string.upper()+'='+str(mag)[:4]+','+comment_string)
    #for other in thing:
        #print(other.replace(['d','h','m','s'], ':'))
print(name_array)

output_array= np.vstack([name_list, ra_list, dec_list, epoch_list, mag_list])
if sort_by_ra:
    sorted_order= np.argsort(output_array[1])
    print(sorted_order)
    print('sorted_order:',output_array[:,sorted_order])
    output_array= output_array[:,sorted_order]
print(output_array)
print(output_array.dtype)
output_array= output_array.T
np.savetxt(output_file,output_array,  delimiter= delimiter, fmt= '%1s')


