"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-09-17


Import wdwarfdate and obtain total age estimates (and cooling age estimates) of all of the white dwarfs from the open clusters. Then save that
output too.




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
from astropy.table import Table, Column, vstack, join, hstack
#import scipy.interpolate as scinterp
import time
import pyvo
import wdwarfdate
import cluster_utilities as cu
#from astroquery.vizier import Vizier

#from astroquery.gaia import Gaia
#from astroquery.xmatch import XMatch
start=time.time()

def elapsed_time():
    
    now= time.time()
    elapsed=now-start
    hours=elapsed // 3600
    minutes = elapsed %3600 //60
    seconds = elapsed %60
    print('Elapsed time:' + str(hours) +' h ' + str(minutes)+ ' m ' + str(seconds) + ' s')
    return

sys.path.append('../')


#input_file='HR24members_GF21maincat_crossmatch_simbadadded.csv'
input_file='HR24members_GF21maincat_crossmatch_simbadadded_mwddadded.csv'

input_table=Table.read(input_file)
hold_table=input_table.copy()
#input_table=input_table[0:5]

#teffs and loggs for hte Gentile Fusillo columns
##teffs=input_table['TeffH']
##teffs_err=input_table['e_TeffH']

##logg=input_table['loggH']
##logg_err=input_table['e_loggH']


#teffs and loggs for the MWDD columns

#find rows that have errobars with their teff and logg



#teffs=input_table['mwdd_teff']
#teffs_err=input_table['mwdd_Dteff']

#logg=input_table['mwdd_logg']
#logg_err=input_table['mwdd_Dlogg']

def get_ages(input_table,param_source='mwdd',model_wd='DA'):
    """
    param_source either "mwdd" for Montreal White Dwarf Database columns or "gf21" for the Gentile Fusillo et al. 2021 columns
    
    """
    if param_source=='mwdd':
        teffs=input_table['mwdd_teff']
        teffs_err=input_table['mwdd_Dteff']

        logg=input_table['mwdd_logg']
        logg_err=input_table['mwdd_Dlogg']
        #mask=input_table['mwdd_Dlogg'].mask
        #unmasked_indices=np.where(~mask)[0]
        #input_table=input_table[unmasked_indices]
        #print('only valid errorbars')
        #input_table.pprint()
        #input_gaiadr3=input_table['GaiaDR3']
    elif param_source=='gf21':
        teffs=input_table['TeffH']
        teffs_err=input_table['e_TeffH']

        logg=input_table['loggH']
        logg_err=input_table['e_loggH']
    
    else:
        print("invalid param_source provided:", param_source)
        print("Don't know how the Teff and Logg columns are formatted therefore")
    #print('get_ages() input_table')
    #for colname in input_table.colnames:
        #print(colname)
    WD=wdwarfdate.WhiteDwarf(teffs,teffs_err, logg, logg_err,model_wd=model_wd,feh='p0.00',vvcrit='0.0',model_ifmr='Cummings_2018_PARSEC',high_perc=84,low_perc=16,datatype='Gyr',display_plots=False)
    print('Evaluating ages for '+model_wd+ ' white dwarfs.\n\n')
    WD.calc_wd_age()
    #print('\n\ncolnames in input_table')
    #for colname in input_table.colnames:
        #print(colname)
    print("post age merge input_table['GaiaDR3']",input_table['GaiaDR3'])
    merged=hstack([input_table['GaiaDR3'],WD.results])
    merged['wdwarfdate_model']=np.full(len(merged),model_wd,dtype='S8')
    merged.pprint()
    return merged


def iterate_spectype_ages(first_table, param_source='mwdd',flexible_spectype=True):
    #model_type_list=['DA','non-DA','non-DA','non-DA','non-DA','DA']
    #spectype_list=['DA','DB','DO','DQ','DZ','DC']
    first_table.add_index('GaiaDR3') #should now be able to find rows using GaiaDR3 IDs direclty without having to do extra work.
    #print('\n\nfirst_table.indices',first_table.indices,'\n\n')
    join_needed=True
    if param_source=='mwdd':
        mask=first_table['mwdd_Dlogg'].mask
        unmasked_indices=np.where(~mask)[0]
        input_table=first_table[unmasked_indices]
        input_gaiadr3=input_table['GaiaDR3']
    else:
        input_table=first_table
    #for colname in input_table.colnames:
        #print(colname)
    model_type_list=['DA','non-DA','non-DA','non-DA','DA']
    spectype_list=['DA','DB','DO','DQ','DC']
    #model_type_list=['non-DA','non-DA','non-DA','DA']
    #spectype_list=['DB','DO','DQ','DC']
    if param_source=='mwdd':
        spectype_colname='mwdd_spectype'
    else:
        spectype_colname='simbad_sp_type'
    for model, spectype in zip(model_type_list,spectype_list):
        print('\n\nStarting analysis for '+spectype+ ' white dwarfs using ' + model+ ' models.')
        if flexible_spectype:
            print('Flexible spectral type allowed.')
        else:
            print('Inflexible spectral type. (Must be exact match).')
        #for row in input_table:
            #print(row['mwdd_spectype']==spectype,row['mwdd_spectype'], spectype)
        subtable=cu.limit_spectype(input_table,spectype=spectype,flexible_spectype=flexible_spectype,spectype_colname=spectype_colname)
        age_table=get_ages(subtable,param_source=param_source,model_wd=model)
        if join_needed:
            try:
                first_table=join(first_table,age_table,keys_left='GaiaDR3',keys_right='GaiaDR3',join_type='outer')
                for colname in first_table.colnames:
                    if (('recno'  in colname) or ('RV' in colname)):
                        pass
                    else:
                        if ('_1' in colname):
                            print('removing _1')
                            first_table.rename_column(colname,colname.replace('_1',''))
                        elif '_2' in colname:
                            print('removing _2')
                            first_table.remove_column(colname)
                        else:
                            pass
                join_needed=False
                for colname in first_table.colnames:
                    print(colname)
                first_table.add_index('GaiaDR3')
                #input_table.rename_column('GaiaDR3_1','GaiaDR3')
                #input_table.remove_column('GaiaDR3_2')
            except ValueError as error:
                print("\n\nValueError:",error)
                print("Probably means we had an empty table for this spectral type:",spectype)
                subtable.pprint()
                if len(subtable)==0:
                    print('Yep, the subtable for this spectral type was empty.')
                else:
                    print('Weird... not sure why that choked then...')
        else:
            print('attempting to slide age_table columns into the existing one.')
            if len(subtable)==0:
                    print('Yep, the subtable for this spectral type was empty.')
            else:
                #first_table=join(first_table,age_table,keys_left='GaiaDR3',keys_right='GaiaDR3',join_type='outer')
                #for colname in first_table.colnames:
                    #if (('recno'  in colname) or ('RV' in colname)):
                        #pass
                    #else:
                        #if ('_1' in colname):
                            #first_table.rename_column(colname,colname.replace('_1',''))
                        #elif '_2' in colname:
                            #first_table.remove_column(colname)
                        #else:
                            #pass
                for colname in age_table.colnames:
                    #age_table.add_index('GaiaDR3')
                    #print(first_table.indices)
                    for num,gaiaid in enumerate(age_table['GaiaDR3']):
                        inds=np.where(first_table['GaiaDR3']==gaiaid)[0]
                        #print(num,inds,gaiaid,first_table[colname][inds],age_table[num][colname])
                        #print(type(num),type(inds),type(gaiaid),type(first_table[colname][inds][0]),type(age_table[num][colname]))
                        for single_index in inds:
                            first_table[colname][single_index]=age_table[num][colname]
                    #first_table.loc[][colname]=age_table[colname]
                    #print("age_table['GaiaDR3']",age_table['GaiaDR3'])
                    #print(first_table.loc[age_table['GaiaDR3']][colname])
        elapsed_time()
    return first_table

#WD=wdwarfdate.WhiteDwarf(teffs,teffs_err, logg, logg_err,model_wd='DA',feh='p0.00',vvcrit='0.0',model_ifmr='Cummings_2018_PARSEC',high_perc=84,low_perc=16,datatype='Gyr',display_plots=False)

#WD.calc_wd_age()

#WD.results

#WD.results.pprint()


##merged=join(input_table,WD.results)

#merged=hstack([input_table,WD.results])

#super_merged=join(hold_table,input_table,keys_left='GaiaDR3', keys_right='GaiaDR3',join_type='left')
#input_table=input_table[:10]
merged=iterate_spectype_ages(input_table,param_source='mwdd', flexible_spectype=True)
output_name=input_file.split('.')[0]+'_wdagesadded.csv'

print('Saving ' + output_name)

merged.write(output_name)

print(output_name + ' saved... or failed to save, but the script executed the save attempt.')






