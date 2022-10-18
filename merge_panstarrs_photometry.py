"""
Created by Ben Kaiser (UNC-Chapel Hill)2022-05-02 from mark_repeats.py

This should be adding columns to the existing Gaia file based on matches to panstarrs 2 that will include the grizY photometry from panstarrs 2 for that object.

"""

from __future__ import print_function
import numpy as np
from astropy.table import Table, QTable, Column

print('\n\n\n\n##################\n\n\n\n\n##############')

#newest_search= '20190405_purple_search_gaia_sc.csv'
#newest_search= '20190516B_retargeted_purple_search_gaia_scbd.csv'
#newest_search= '20190917_alkaliWD_attempt2_gaia_scbd.csv'
#check_searches=['20190107_chris_merge_gaia.csv',
                #'expanded_purple_search_gmaglimit_gaia_sc.csv',
                #'exc1_8_2_2_purple_search_gmaglimit_gaia_sc.csv',
                #'sdssj1330_similar_subset_gaia_sc.csv',
                #'20190405_purple_search_gaia_sc.csv',
                #'20190829_alkaliWD_targeted_gaia_scbd.csv',
                #'20190516B_retargeted_purple_search_gaia_scbd.csv'] #old filenames that should be searched for duplicates


newest_search='20190516B_retargeted_purple_search_gaia_scbd_20220503_update.csv'

#check_searches=['20200819_retargeted_purple_search_wpanstarrs_REALLYjust_ra_dec_csv_5_2_2022_panstarrs2_params_properRAsort_goodmags.csv']

#check_searches=['20200819_retargeted_purple_search_wpanstarrs_REALLYjust_ra_dec_csv_5_3_2022_panstarrs1_params_properRAsort_goodmags.csv']
check_searches=['20200819_retargeted_purple_search_wpanstarrs_REALLYjust_ra_dec_csv_5_3_2022_panstarrs1stack_params_properRAsort_goodmags.csv']

comp_name= 'source_id' #column name that should be used to check if it's a repeat Probably should be 'source_id' since that's it's Gaia-assigned name. A long string of numbers unique to each object.

color_list=['g','r','i','z','y']
#color_string_base='MeanPSFMag'
color_string_base='PSFMag'


#create the new column to be occupied with indications of repeat. Initialized with False for every entry.
newest_table= Table.read(newest_search)



for color in color_list:
    
    try:
        print('Test',color , 'column: ', newest_table[color+color_string_base][0])
    except KeyError:
        print(color, " column doesn't exist, so we're making it.")
        new_length = len(newest_table['dist'])
        false_array = np.empty((new_length,), dtype=float)
        false_array[:]=np.nan
        #false_column = Column(false_array,name='ps2_'+color+'_mean_mag',dtype=float)
        #newest_table.add_column(false_column)
        #false_column = Column(false_array, name='ps2_'+color+'_mean_mag_error',dtype=float)
        #newest_table.add_column(false_column)
        false_column = Column(false_array,name='ps1_'+color+'_mean_mag',dtype=float)
        newest_table.add_column(false_column)
        false_column = Column(false_array, name='ps1_'+color+'_mean_mag_error',dtype=float)
        newest_table.add_column(false_column)

#try:
    #print('Test prev_file column: ', newest_table['prev_file'][0])
#except KeyError:
    #print("prev_file column doesn't exist, so we're making it.")
    #new_length = len(newest_table['dist'])
    #file_array = np.full((new_length,),'',  dtype='S64')
    #file_column = Column(file_array, name='prev_file',dtype='S64')
    #newest_table.add_column(file_column)

#newest_table['name']=newest_table['name'].astype(np.str_)
print(newest_table['name'].dtype)
newest_table.pprint()

def mark_repeats(old_table, new_table, old_file):
    for row in new_table:
        #print(comp_entry)
        #old_comp= np.copy(comp_entry)
        #comp_entry=comp_entry.astype(np.int_)
        #print(comp_entry)
        #print(old_comp-comp_entry)
        #matched_index= np.where(new_table[comp_name].astype(np.int_)==comp_entry.astype(np.int_))
        
        #matched_index= np.where(np.abs(row['ra']-old_table['_ra_'])+np.abs(row['dec']-old_table['_dec_']) < 0.0000001)
        #matched_index= np.where(np.abs(row['ra']-old_table['_ra_'])+np.abs(row['dec']-old_table['_dec_']) < 0.01)
        #print('matched_index', matched_index)
        #print('matching dstArcSec of ps2', old_table['dstArcSec'][matched_index])
        #min_dist_order=np.argsort(old_table['dstArcSec'][matched_index])
        #print('\n\n++++',min_dist_order,'++++\n\n')
        
        separations=np.sqrt((row['ra']-old_table['MatchRA'])**2+(row['dec']-old_table['MatchDEC'])**2)
        min_index=np.argmin(separations)
        print('min separation:', separations[min_index], row['ra'], old_table['MatchRA'][min_index], row['dec'], old_table['MatchDEC'][min_index])
        if separations[min_index] > 10.:
            min_index=np.nan
            print('Big separation\n----')
        else:
            pass
        try:
            #print('comp_entry.dtype', comp_entry.dtype)
            #print('new_table[comp_name].dtype', new_table[comp_name].dtype)
            #print("new_table['repeat'][matched_index].dtype", new_table['repeat'][matched_index].dtype)
            #print("new_table['name'][matched_index].dtype", new_table['name'][matched_index].dtype)
            #new_table['repeat'][matched_index]=True
            #new_table['name'][matched_index]=name
            #print('old_file', old_file)
            #new_table['prev_file'][matched_index]=old_file
            #old_table_row=old_table[matched_index][min_dist_order][0]
            old_table_row=old_table[min_index]
            
            print('old_table_row')
            print(old_table_row['dstArcSec'])
            #print(old_table_row)
            #old_table_row.pprint()
            for color in color_list:
                #row['ps2_'+color+'_mean_mag']=old_table[color+color_string_base][matched_index]
                #row['ps2_'+color+'_mean_mag_error']=old_table[color+color_string_base+'Err'][matched_index]
                #row['ps2_'+color+'_mean_mag']=old_table_row[color+color_string_base]
                #row['ps2_'+color+'_mean_mag_error']=old_table_row[color+color_string_base+'Err']
                row['ps1_'+color+'_mean_mag']=old_table_row[color+color_string_base]
                row['ps1_'+color+'_mean_mag_error']=old_table_row[color+color_string_base+'Err']
        except ValueError as error:
            print('ValueError',error)
            pass
        except IndexError as error:
            print('IndexError', error)
    return new_table


for old_file in check_searches:
    old_table = Table.read(old_file)
    done_table = mark_repeats(old_table, newest_table, old_file)
    done_table.pprint()
    

#done_table.write(newest_search, format= 'ascii.csv', overwrite=True)

#new_filename=newest_search.split('.')[0]+'_ps2_phot_added.csv'

new_filename=newest_search.split('.')[0]+'_ps1stack_phot_added.csv'

done_table.write(new_filename,format='ascii.csv')





