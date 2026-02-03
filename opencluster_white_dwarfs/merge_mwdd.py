"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-09-18

This should ingest the existing master table of open cluster white dwarfs and an ID-based cross match I performed on them in the MWDD portal and then join them together to a new master table.

I'm also going to have it clean up the existing columns by removing all columns with _number indicating a repeat of information. That should at least partially fix the problem.

I'm also going to add the new columns with "mwdd_" stuck onto the front of them so it's obvious where they came from.




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
from astropy.table import Table, Column, vstack, join, unique
#import scipy.interpolate as scinterp
import time
import pyvo


sys.path.append('../')

existing_file='HR24members_GF21maincat_crossmatch_simbadadded.csv' #I'm not using the one that already has the wdwarfdate done because I'm planning to use these Teff and logg parameters to get it.
mwdd_file='mwdd_openclustermatches.csv'


existing_table=Table.read(existing_file)
mwdd_table=Table.read(mwdd_file)

for colname in existing_table.colnames:
    print('colname',colname)
    if (('recno' in colname) or ('RV' in colname)):
        pass
    else:
        if ('_1' in colname):
            print('colname with _1',colname)
            print('new name:',colname.replace('_1',''))
            print('comparison of values for 1 vs 2 column', existing_table[colname]==existing_table[colname.replace('_1','_2')])
            print('side-by-side comparison:')
            existing_table[colname,colname.replace('_1','_2')].pprint()
                                                                                    
            try:
                print('existing col?:',existing_table[colname.replace('_1','')])
            except KeyError:
                print("that column doesn't already exist, yay!\nSure hope the value isn't the same for the second one... maybe I should just check that...")
            try:
                existing_table.rename_column(colname,colname.replace('_1',''))
            except KeyError as error:
                print("KeyError:",error)
                print("Column already exists, so we're going to discard the _1 version")
                existing_table.remove_column(colname)
        elif('_2' in colname):
            existing_table.remove_column(colname)
    else:
        print("recno for some reason isn't always the same... so we're not going to remove it or change its name.")
        
        
mwdd_table.rename_column('request','gaiaedr3') #For some reason the "Gaia DR3" field on the MWDD doesn't populate for every object for which it actually possesses a Gaia DR3 identifier... so this is based on the number I uploaded that it still matched to one of the various names for each object (in this case the DR3 identifier almost certainly)
for colname in mwdd_table.colnames:
    mwdd_table.rename_column(colname,'mwdd_'+colname)


mwdd_table['mwdd_gaiaedr3'].pprint()

merged=join(existing_table,mwdd_table, keys_left='GaiaDR3',keys_right='mwdd_gaiaedr3', join_type='left')
merged.pprint()
print(len(merged),len(mwdd_table),len(existing_table))

unique_nums, unique_inds=np.unique(merged['mwdd_gaiaedr3'],return_index=True)
for index in unique_inds:
    inds=np.where(merged['mwdd_gaiaedr3']==merged['mwdd_gaiaedr3'][index])
    print(len(inds[0]),inds)
    if len(inds[0])>1:
        merged[inds].pprint()
        sum_matches=merged[inds[0][0]]==merged[inds[0][1]]
        print(sum_matches)
        print(np.sum(sum_matches))
        print(merged['simbad_sp_type'][inds])
        print(merged['simbad_name'][inds])
        
merged_unique=unique(merged,keys='GaiaDR3',keep='first') #remove duplicates. I determined the duplicates were sdBs so it doesn't matter anyway. I guess they have more than one entry in the MWDD.
print(merged_unique)


output_name=existing_file.split('.')[0]+'_mwddadded.csv'
#merged.write(output_name)
merged_unique.write(output_name)
