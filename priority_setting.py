"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-05-28

This script should take one of my csv versions of gaia results and append a new column to the astropy table
that sets the priority level of it for observing based on numerical criteria, so I don't have to individually go through 
them one by one assigning priority values.


"""
from __future__ import print_function
import numpy as np

import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column
import matplotlib.pyplot as plt
import time
start = time.time()

input_filename= '20190516B_retargeted_purple_search_gaia_scbd.csv'

output_filename= input_filename #I want to overwrite the input file, so I don't have proliferating files

input_table = Table.read(input_filename)


base_priority = 2

##################################
table_length= len(input_table['ra']) #establish the length of the table

priority_column= Column(np.ones(table_length), name="priority",dtype=int)
priority_column=priority_column*base_priority
try:
    print(input_table['priority'])
    print("priority column already exists, so we'll be overwriting those values")
    input_table['priority']=priority_column
except KeyError as error:
    print("KeyError:", error)
    print("priority column does not yet exist, so we'll create it.")
    input_table.add_column(priority_column)

#############################


def set_pwd_priority(table, rescale_val=4):
    """
    adjust the 'priority' column of the table based on its pwd value
    
    """
    additional_points = np.int_(0.9999+ table['pwd']*rescale_val) #adding 0.999 should make it so non-zero Pwd values now gain a point
    table['priority']=table['priority']+additional_points
    return table

def set_g_rp_priority(table, cutoff=1.00):
    """
    adjust the 'priority column of the table based on its g_rp value
    
    if g_rp < cutoff, then priority= priority +1
    """
    add_indices= np.where(table['g_rp']<cutoff)
    add_array = np.zeros(table_length)
    add_array[add_indices]= 1
    table['priority']=table['priority']+add_array
    return table


##################################

input_table = set_pwd_priority(input_table)
input_table = set_g_rp_priority(input_table)

input_table['priority']=np.int_(input_table['priority']) #making sure this column is an integer value, hopefully

input_table.pprint()

print('saving table to ', output_filename)
input_table.write(output_filename, format='ascii.csv', overwrite=True)
print('table saved')

plt.scatter( input_table['g_rp'], input_table['pwd'],c=input_table['g_rp'])
plt.ylabel('Pwd')
plt.xlabel('g_rp')
plt.show()

plt.scatter(input_table['pwd'], input_table['priority'], c=input_table['g_rp'])
plt.xlabel('Pwd')
plt.ylabel('Priority')
plt.show()



plt.scatter( input_table['g_rp'], input_table['pwd'],c=np.array(input_table['priority']))
plt.ylabel('Pwd')
plt.xlabel('g_rp')
plt.show()
