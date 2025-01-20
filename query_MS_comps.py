"""
Created by Ben Kaiser (UNC-Chapel Hill) 2024-01-30

This should take the dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust.fits file and 
compare that with the white dwarf coordinates (or eDR3 ID or DR2 ID) from 
dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf_Vincent2023bClass_20240125update.csv. It should then maybe query the MS DR3 data out of Gaia 
perhaps and maybe append it to the same table?

Actually I should probably make a separate output, but to that I should include the WD RA 
and Dec  and Gmag so I can recover my white dwarfs in this other table. 

This file will borrow heavily from queryDR3_from_DR2file.py

"""


import numpy as np
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column
import matplotlib.pyplot as plt

import time
import sys
start = time.time()

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" 

input_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust.fits'
wd_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf_Vincent2023bClass_20240125update.csv'

output_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_MSmatchesIDs.csv'
    
if Gaia.MAIN_GAIA_TABLE=='gaiadr3.gaia_source':
    gaia_string='gaiaDR3'
elif Gaia.MAIN_GAIA_TABLE=='gaiaedr3.gaia_source':
    gaia_string='gaiaeDR3'
elif Gaia.MAIN_GAIA_TABLE=='gaiadr2.gaia_source':
    gaia_string='gaiaDR2'
else:
    pass


credentials_file= 'Gaia_credentials.txt'
print("Logging in.")
Gaia.login(credentials_file= credentials_file)
print("Log in succesful.")


coord_list= []


all_table=Table.read(input_file)

wd_table=Table.read(wd_file,delimiter=',')


all_ms_ra=np.array([])
all_ms_dec=np.array([])
all_wd_ra=np.array([])
all_wd_dec=np.array([])


table_initialized=False
#for row in output_fulltable:
    #wd_num=row['wd_comp']
    #if table_initialized:
        #sub_row=row[['ra'+str(wd_num),'dec'+str(wd_num)]]
        #wd_output_table.add_row(sub_row)
    #else:
        #sub_row=row[['ra'+str(wd_num),'dec'+str(wd_num)]]
        #wd_output_table=Table(names=new_colnames,rows=sub_row)
        #table_initialized=True
count=0
new_colnames=['ra','dec','phot_g_mean_mag','pairdistance','sep_AU','R_chance_align','ra_wd','dec_wd','phot_g_mean_mag_wd']
new_colnames=['source_id','ra','dec','phot_g_mean_mag','pairdistance','sep_AU','R_chance_align','source_id_wd','ra_wd','dec_wd','phot_g_mean_mag_wd']

for row in all_table:
    wd_ra=row['ra'+str(row['wd_comp'])]
    wd_dec=row['dec'+str(row['wd_comp'])]
    wd_num=row['wd_comp']
    if row['wd_comp']==1:
        ms_comp=2
    elif row['wd_comp']==2:
        ms_comp=1
    else:
        print('\n\nsomehow WD comp is neither 1 nor 2?\n\n')

    matched_index= np.where(np.abs(wd_table['ra']-wd_ra)+np.abs(wd_table['dec']-wd_dec) < 0.0000001)
    if matched_index[0].shape[0]>0:
        wd_row=wd_table[matched_index[0]]
        print('matched_index', matched_index)
        sub_row=row[['source_id'+str(ms_comp),
            'ra'+str(ms_comp),'dec'+str(ms_comp),'phot_g_mean_mag'+str(ms_comp),
                'pairdistance','sep_AU','R_chance_align',
                'source_id'+str(wd_num),
                'ra'+str(wd_num),'dec'+str(wd_num),'phot_g_mean_mag'+str(wd_num)]]
        if table_initialized:
            #sub_row=row[['ra'+str(ms_comp),'dec'+str(ms_comp),'phot_g_mean_mag'+str(ms_comp),
                         #'pairdistance','sep_AU','R_chance_align',
                         #'ra'+str(wd_num),'dec'+str(wd_num),'phot_g_mean_mag'+str(wd_num)]]
             #sub_row=row[['source_id'+str(ms_comp)
                 #,'ra'+str(ms_comp),'dec'+str(ms_comp),'phot_g_mean_mag'+str(ms_comp),
                         #'pairdistance','sep_AU','R_chance_align',
                         #'source_id'+str(wd_num)
                         #'ra'+str(wd_num),'dec'+str(wd_num),'phot_g_mean_mag'+str(wd_num)]]
            output_table.add_row(sub_row)
        else:
            #sub_row=row[['ra'+str(ms_comp),'dec'+str(ms_comp),'phot_g_mean_mag'+str(ms_comp),
                         #'pairdistance','sep_AU','R_chance_align',
                         #'ra'+str(wd_num),'dec'+str(wd_num),'phot_g_mean_mag'+str(wd_num)]]
            output_table=Table(names=new_colnames,rows=sub_row)
            table_initialized=True
        print('\ndif RA before',output_table[count]['ra_wd']-wd_row['ra'])
        output_table[count]['ra_wd']=wd_row['ra']
        output_table[count]['dec_wd']=wd_row['dec']
        output_table[count]['phot_g_mean_mag_wd']=wd_row['phot_g_mean_mag']
        print('dif RA after',output_table[count]['ra_wd']-wd_row['ra'],'\n')
        
        count+=1
    
    else:
        pass

print(count)
print(len(wd_table))

output_table.pprint()

output_table.write(output_file,format='ascii.csv')




