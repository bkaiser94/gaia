"""
Created by Ben Kaiser (UNC-Chapel Hill) 2024-03-15.

This script should take the combined SDSS-Gaia table of the mistaken K dwarfs and actually just directly retrieve all of the SDSS spectra from it. 
I might have it also rename the spectra appropriately, but I'm not sure.



"""
from __future__ import print_function

from astroquery.sdss import SDSS

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.table import join as ATjoin
from astropy.table import hstack as AThstack
from astropy.table import vstack as ATvstack
from astropy import units as u
from astropy import constants as const
from astropy import coordinates as coords
from astropy.io import fits
import sys
import os

input_file='WDs_mistaken_for_Kstars_nofehlimit_by_SDSS_gaiaDR3.csv'

output_dir='/Users/BenKaiser/Desktop/SDSS_speclib/Kstars_mistaken_for_WDs/'

input_table=Table.read(input_file)
default_name_base='SDSSJ'

#sdss_data_release=14
#sdss_data_release=17
#sdss_data_release=13
#sdss_data_release=11 #doesn't have astroquery field info
sdss_data_release=12 #"penultimate SDSS-III data release according to its website"
#sdss_data_release=7 #doesn't have astroquery field info

#positions=coords.SkyCoord(ra=input_table['sdss_RA'],dec= input_table['sdss_DEC'],frame='icrs', unit=(u.deg, u.deg))

#print(len(positions))

#xid=SDSS.query_region(positions, radius='1 arcsec')
#print("xid complete")
#print(xid)
#print(len(xid))
def make_name(row,filename=True, name_base=default_name_base):
    pos=coords.SkyCoord(ra=row['ra'], dec=row['dec'], unit=(u.deg, u.deg))
    string_coords=pos.to_string(style='hmsdms')
    print(type(string_coords))
    print(string_coords)
    replace_chars= ['d','h','m']
    for char in replace_chars:
        string_coords= string_coords.replace(char, '')
    string_coords=string_coords.replace('s', '')
    split_string=string_coords.split(' ')
    ra= split_string[0]#limiting precision of the decimal
    dec= split_string[1] #limiting precision of the decimal, need the additional index for sign
    small_ra= ra[:4]
    small_dec=dec[:5] #need additional index for + or -
    print('small_dec', small_dec)
    if filename:
        small_dec=small_dec.replace('-','m')
        small_dec=small_dec.replace('+','p')
        print('small_dec', small_dec)
    else:
        pass
    name=name_base+small_ra+small_dec
    return name


count=0
failed_count=0
existing_count=0
new_retrievals=0
for row in input_table:
    print('\n===========')
    print('Retrieving spec ', count)
    print('Plate:', row['plate'], 'MJD:', row['mjd'], 'FiberID:', row['fiberID'])
    filename_core=make_name(row, filename=True)
    filename=filename_core+'_sdss_spec.fits'
    filename=output_dir+filename
    if os.path.exists(filename):
        print(filename,'exists, so skipping the shenanigans')
        existing_count+=1
    else:
        try:
            spec=SDSS.get_spectra(plate=row['plate'], mjd=row['mjd'], fiberID=row['fiberID'],data_release=sdss_data_release)
            #print(spec)
        
            print('saving to', filename)
            print('elements of spec')
            #for element in spec:
                #print(element, type(element))
                
            hdulist=fits.HDUList(spec[0])
            hdulist.writeto(filename)
            print('file saved to', filename)
            new_retrievals+=1
        except OSError as error:
            print('OSError:', error)
            print('So skipping it and continuing')
            existing_count+=1
        except TypeError as error:
            print('TypeError:', error)
            print("I guess that means this object didn't have a valid spectrum returned. This was row "+str(count)+' as a reminder.')
            print('Plate:', row['plate'], 'MJD:', row['mjd'], 'FiberID:', row['fiberID'])
            print(filename_core)
            print('Failed to save', filename)
            failed_count+=1
    
    
    
    
    count+=1
    
print('Number of spectra that could not be retrieved:', failed_count)
print('Number of newly retrieved and saved spectra:', new_retrievals)
print('Number of spectra that had already been downloaded:', existing_count)

print('Total number of specta dealt with then (sum of those numbers):', failed_count+new_retrievals+existing_count)
print("length of input table (should be same number as above):", len(input_table))
