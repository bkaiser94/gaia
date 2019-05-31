"""
Created by Ben Kaiser 2019-05-10 (UNC-Chapel Hill)


Take multiple target list inputs (already formatted for SOAR more or less, meaning with or without target 
numbers), then put them all together in order as provided (so don't reorder by RA or anything), and then
put numbers leading the names beginning at a number specified.


"""
from __future__ import print_function
import numpy as np


starting_num=700

#output_filename='WDdMcands_target_list.txt'
#output_filename='coolWDcands_target_list.txt'
output_filename='20190516B_retargeted_purple_search_gaia_scbd_targlist.txt'

list_of_names=['20190516B_retargeted_purple_search_gaia_scbd_targlist.txt']
#list_of_names=['sdssj1330_similar_subset_gaia_sc_targlist.txt']

#list_of_names =['20190107_chris_merge_gaia_targlist_wcomps.txt',
                #'expanded_purple_search_gmaglimit_wcomps_targlist.txt',
                #'exc1_8_2_2_purple_search_gmaglimit_wcomps_targlist.txt',
                #'20190405_purple_search_gaia_unique_wcomps_targlist.txt'
                #]
#list_of_names =[
                #'expanded_purple_search_gmaglimit_wcomps_targlist.txt',
                #'exc1_8_2_2_purple_search_gmaglimit_wcomps_targlist.txt',
                #'20190405_purple_search_gaia_unique_wcomps_targlist.txt'
                #]

list_of_lists=[]

for name in list_of_names:
    print('name',name)
    input_array = np.genfromtxt(name, delimiter='\t',dtype=np.str).T
    #input_array = np.genfromtxt(name, delimiter='\t')
    #print('input_array.shape', input_array.shape)
    print(input_array)
    list_of_lists.append(input_array)
    
print(list_of_lists)

array_of_lists= np.hstack(list_of_lists)
print(array_of_lists.shape)
print(array_of_lists[0])

number_range= np.arange(starting_num, starting_num+array_of_lists.shape[1], 1)
number_range=number_range.astype(np.str)
print(number_range)
number_range=number_range.tolist()
first_col=array_of_lists[0].tolist()
print(first_col)
#np.str_(number_range)+ array_of_lists[0]
#first_col= '_'.join(zip(number_range, array_of_lists[0]))
first_col= map('_'.join,zip(number_range, first_col))
#first_col=np.array([x1 + '_'+ x2 for x1,x2 in zip(number_range,array_of_lists[0])])
#print(first_col.shape)
print(first_col)

#array_of_lists[0]=first_col
array_of_lists=np.vstack([first_col,array_of_lists[1:]])

print(array_of_lists)

np.savetxt(output_filename, array_of_lists.T, delimiter='\t', fmt='%1s')
