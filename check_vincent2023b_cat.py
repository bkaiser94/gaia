"""
Created by Ben Kaiser (UNC-Chapel Hill) 2023-09-12

This script should take whatever input set of Gaia objects (probably DR3 generally if not always) and then check Vincent et al. 2023b's catalogue of
 Machine-learning classified XP spectra for high probability white dwarfs to see if they provide a classfication for that object.

I'll probably just append the columns to the original table or something. I don't know I haven't thought that far yet.



"""


import numpy as np
#from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column, join
import time
import sys
start = time.time()



vincent_file='GSPCWD_catalogue.csv'

input_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf.csv'

output_file=input_file.split('.')[0]+'_Vincent2023bClass.csv'


########################


input_table=Table.read(input_file)

vincent_table=Table.read(vincent_file)

joined_table=join(input_table,vincent_table,keys='source_id',join_type='left')
joined_table.pprint()

#joined_table.write(output_file,overwrite=True)
sys.exit()

"""
Everything past this part is actually superfluous to incorporating the data from the Vincent table. Turns out the join function of astropy.table has
been able to handle this all along. I suppose the stuff down below could also used to get the count of different types of objects, but that really doesn't 
seem necessary.


"""

############################################
input_table.pprint()


vincent_table.pprint()

plist=['P_DA','P_DB','P_DC','P_DO','P_DQ','P_DZ']

print(vincent_table.dtype)

vincent_colnames=vincent_table.colnames[1:] #excludes the source_id as it is an unnecessary repeat
print('len input_table', len(input_table))
print(len(vincent_colnames))
print(vincent_colnames)
print(type(vincent_colnames))
print('zeros array',np.zeros((len(vincent_colnames),len(input_table)),dtype='S168'))
print(np.zeros((len(vincent_colnames),len(input_table)),dtype='S168').shape)
input_table.add_columns(np.zeros((len(vincent_colnames),len(input_table)),dtype='S168'),names=vincent_colnames)
for name in vincent_colnames:
    print(input_table[name].dtype, vincent_table[name].dtype)
    input_table[name].dtype=vincent_table[name].dtype
    print('len name',len(input_table[name]))
count_dict={
    'DA':0,
    'DB':0,
    'DC':0,
    'DO':0,
    'DQ':0,
    'DZ':0
    }
#DA_count=0
#DB_count=0
#DC_count=0
#DO_count=0
#DQ_count=0
#DZ_count=0
spec_list=[]

for row in input_table:
    matched_index=np.where(vincent_table['source_id']==row['source_id'])
    print('matched_index',matched_index)
    try:
        matched_row=vincent_table[matched_index]
        #row[vincent_colnames]=matched_row[vincent_colnames]
        for name in vincent_colnames:
            #print('row['+name+']', row[name])
            #print('matched_row['+name+'][0]',matched_row[name][0])
            
            row[name]=matched_row[name][0]
            #print(type(row[name])==type(matched_row[name][0]))
        #print(len(matched_row[vincent_colnames]))
        #print('row[vincent_colnames]',row[vincent_colnames])
        #print('matched_row[vincent_colnames]',matched_row[vincent_colnames])
        spectype=str(matched_row['spectype'][0])
        for name in count_dict:
            if name in spectype:
                count_dict[name]+=1
        
        print('spectral type',matched_row['spectype'])
        if 'DZ' in spectype:
            for prob in plist:
                print(matched_row[prob])
        spec_list.append(matched_row['spectype'])
    except IndexError as error:
        print("IndexError so it didn't try to do the addition of values")
    print('\n\n=============')
    
#input_table[vincent_colnames].dtype=vincent_table[vincent_colnames].dtype
print(count_dict)
    
    
input_table.pprint()

for name in vincent_colnames:
    try:
        bad_vals=np.where(np.abs(input_table[name])<1e-20)
        input_table[name][bad_vals]=0
    except np.core._exceptions.UFuncTypeError:
        pass
    

    
input_table.pprint()


for name in input_table.colnames:
    print(name,len(input_table[name]))
    
count=0
for row in input_table:
    count+=1
    print(row)    

print('total rows:',count)

#input_table.write(output_file)










    
