"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-04-05

Creates a new column, if it doesn't already exist, called "repeat" that is used to indicate if this target has already 
been identified previously. It should then replace the name of the object with the other one that wasin the old
version of the list.

"""

from __future__ import print_function
import numpy as np
from astropy.table import Table, QTable, Column

print('\n\n\n\n##################\n\n\n\n\n##############')

#newest_search= '20190405_purple_search_gaia_sc.csv'
#newest_search= '20190516B_retargeted_purple_search_gaia_scbd.csv'
newest_search= '20190917_alkaliWD_attempt2_gaia_scbd.csv'
check_searches=['20190107_chris_merge_gaia.csv',
                'expanded_purple_search_gmaglimit_gaia_sc.csv',
                'exc1_8_2_2_purple_search_gmaglimit_gaia_sc.csv',
                'sdssj1330_similar_subset_gaia_sc.csv',
                '20190405_purple_search_gaia_sc.csv',
                '20190829_alkaliWD_targeted_gaia_scbd.csv',
                '20190516B_retargeted_purple_search_gaia_scbd.csv'] #old filenames that should be searched for duplicates

comp_name= 'source_id' #column name that should be used to check if it's a repeat Probably should be 'source_id' since that's it's Gaia-assigned name. A long string of numbers unique to each object.


#create the new column to be occupied with indications of repeat. Initialized with False for every entry.
newest_table= Table.read(newest_search)


try:
    print('Test repeat column: ', newest_table['repeat'][0])
except KeyError:
    print("Repeat column doesn't exist, so we're making it.")
    new_length = len(newest_table['dist'])
    false_array = np.full((new_length,), False, dtype=bool)
    false_column = Column(false_array, name='repeat',dtype=bool)
    newest_table.add_column(false_column)

try:
    print('Test prev_file column: ', newest_table['prev_file'][0])
except KeyError:
    print("prev_file column doesn't exist, so we're making it.")
    new_length = len(newest_table['dist'])
    file_array = np.full((new_length,),'',  dtype='S64')
    file_column = Column(file_array, name='prev_file',dtype='S64')
    newest_table.add_column(file_column)

newest_table['name']=newest_table['name'].astype(np.str_)
print(newest_table['name'].dtype)
newest_table.pprint()

def mark_repeats(old_table, new_table, old_file):
    for ra, dec, name in zip(old_table['ra'], old_table['dec'], old_table['name']):
        #print(comp_entry)
        #old_comp= np.copy(comp_entry)
        #comp_entry=comp_entry.astype(np.int_)
        #print(comp_entry)
        #print(old_comp-comp_entry)
        #matched_index= np.where(new_table[comp_name].astype(np.int_)==comp_entry.astype(np.int_))
        matched_index= np.where(np.abs(new_table['ra']-ra)+np.abs(new_table['dec']-dec) < 0.0000001)
        print('matched_index', matched_index)
        try:
            #print('comp_entry.dtype', comp_entry.dtype)
            #print('new_table[comp_name].dtype', new_table[comp_name].dtype)
            #print("new_table['repeat'][matched_index].dtype", new_table['repeat'][matched_index].dtype)
            #print("new_table['name'][matched_index].dtype", new_table['name'][matched_index].dtype)
            new_table['repeat'][matched_index]=True
            new_table['name'][matched_index]=name
            print('old_file', old_file)
            new_table['prev_file'][matched_index]=old_file
        except ValueError:
            pass
    return new_table


for old_file in check_searches:
    old_table = Table.read(old_file)
    done_table = mark_repeats(old_table, newest_table, old_file)
    done_table.pprint()
    

done_table.write(newest_search, format= 'ascii.csv', overwrite=True)

