"""
Created by Ben Kaiser (UNC-Chapel Hill) 2024-03-18


Create a new column (if it doesn't already exist) to be populated with the names of objects in the file.




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
from astropy.table import Table, Column
#import scipy.interpolate as scinterp
import time



input_filename='WDs_mistaken_for_Kstars_nofehlimit_by_SDSS_gaiaDR3.csv'

output_filename='WDs_mistaken_for_Kstars_nofehlimit_by_SDSS_gaiaDR3_wnames.csv'


object_name_base='SDSS J'


input_table=Table.read(input_filename)
new_table=input_table.copy()

try:
    name_array= np.array(new_table['name'])
except KeyError as error:
    print(error)
    print('No names column, so just making all Gaia names')
    name_array=np.full(new_table['ra'].shape, '')
    new_table.add_column(name_array,name= 'name')


def make_names(this_table):
    this_table['name']=this_table['name'].astype('S32')
    count=0
    for row in this_table:
        input_coords = coord.SkyCoord(ra = row['ra'], dec =row['dec'], unit = (u.deg, u.deg), frame = 'icrs')
        string_coords= input_coords.to_string(style='hmsdms')
        replace_chars= ['d','h','m']
        ra_list=[]
        dec_list=[]
        epoch_list=[]
        mag_list=[]
        name_list=[]
        thing=string_coords
        #print(name)
        #print(thing)
        for char in replace_chars:
            thing= thing.replace(char, ':')
        thing=thing.replace('s', '')
        split_string=thing.split(' ')
        ra= split_string[0][:11] #limiting precision of the decimal
        dec= split_string[1][:12] #limiting precision of the decimal, need the additional index for sign
        small_ra= ra.replace(':', '')[:4]
        small_dec=dec.replace(':','')[:5] #need additional index for + or -
        if 'GaiaJ' in row['name']:
            print("Gaia detected", row['name'])
            row['name']=row['name'].replace('GaiaJ',object_name_base)
            #row['name']=row['name'].replace('-','–')
        elif 'SDSSJ' in row['name']:
            print("SDSS detected", row['name'])
            row['name']=row['name'].replace('SDSSJ','SDSS J')
            #row['name']=row['name'].replace('-','–')
        elif 'WISEA0' in row['name']:
            print("WISEA detected", row['name'])
            row['name']=row['name'].replace('WISEA','WISEA J')
            #row['name']=row['name'].replace('-','–')
        elif 'WISEAJ' in row['name']:
            print("WISEAJ detected", row['name'])
            row['name']=row['name'].replace('WISEAJ','WISEA J')
            #row['name']=row['name'].replace('-','–')
        elif 'PSRJ' in row['name']:
            print("PSRJ detected", row['name'])
            row['name']=row['name'].replace('PSRJ','PSR J')
            #row['name']=row['name'].replace('-','–')
        elif 'LPSMJ' in row['name']:
            print("LPSMJ detected", row['name'])
            row['name']=row['name'].replace('LPSMJ','LSPM J')
            #row['name']=row['name'].replace('-','–')
        elif 'ULASJ' in row['name']:
            print("ULASJ detected", row['name'])
            row['name']=row['name'].replace('ULASJ','ULAS J')
            #row['name']=row['name'].replace('-','–')
        elif 'LEHPM' in row['name']:
            print("LEHPM detected", row['name'])
            row['name']=row['name'].replace('LEHPM','LEHPM ')
        elif 'WDJ' in row['name']:
            row['name']=row['name'].replace('WDJ', 'WD J')
        elif 'LP' in row['name']:
            row['name']=row['name'].replace('LP', 'LP ')
        elif 'ULASJ' in row['name']:
            row['name']=row['name'].replace('ULASJ', 'ULAS J')
        elif row['name'] == "":
            row['name']=object_name_base+small_ra+small_dec

        #row['name']=row['name'].replace('-','–')
        #if name=='.':
            #name=object_name_base+small_ra+small_dec
            #print('new name:', name)
        #elif name=='none':
            #name=object_name_base+small_ra+small_dec
            #print('new name:', name)
        #elif name == '0':
            #name=object_name_base+small_ra+small_dec
            #print('new name:', name)
        #print(thing)
        count+=1
    return this_table

output_table=make_names(new_table)

output_table.write(output_filename, format='ascii.csv')








