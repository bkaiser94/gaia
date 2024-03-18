"""
Created by Ben Kaiser (UNC-Chapel Hill) 2023-10-05

This should take a down-selection from el-Badry's wide binary catalogue and compare it to a previous down-selection from the same catalogue.
The previous selection should have stricter cuts than the new one, so this script will pare down the new catalogue to only those objects that are unique
compared to the original down selection. That way I am effectively searching the parameter space outside the original catalogue without having to 
figure out the boolean version of this query. It would have been bad since they are binaries and either binary component could be the problematic one.



"""


from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.table import join as ATjoin
from astropy.table import hstack as AThstack
from astropy.table import vstack as ATvstack
from astropy import units as u
from astropy import constants as const
import sys




tight_cut_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust.fits'
new_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_loosecuts.fits'
#output_filename='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_loosecutsonly.fits'
output_filename='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_loosecutsonly_justWDradec.csv'


tight_cut_table=Table.read(tight_cut_file)
new_table=Table.read(new_file)



repeat_array=np.zeros(len(new_table))

tally=0
i=0
for row in new_table:
    if row['source_id1'] in tight_cut_table['source_id1']:
        print('source_id1 found in both tables')
        repeat_array[i]=1
        tally+=1
    else:
        pass
    i+=1

print('summed repeats (total repeats):' ,np.sum(repeat_array))
print('tight_cut_table length:',len(tight_cut_table))
print('difference in table lengths:', len(new_table)-len(tight_cut_table))
print('total tally:',tally)

new_indices=np.where(np.int_(repeat_array)==0)
output_fulltable=new_table[new_indices]


output_fulltable.pprint()

#output_table.write(output_fullfilename,overwrite=True)



#sys.exit()
new_colnames=['ra','dec']
print('new_colnames',new_colnames)



table_initialized=False
for row in output_fulltable:
    wd_num=row['wd_comp']
    if table_initialized:
        sub_row=row[['ra'+str(wd_num),'dec'+str(wd_num)]]
        wd_output_table.add_row(sub_row)
    else:
        sub_row=row[['ra'+str(wd_num),'dec'+str(wd_num)]]
        wd_output_table=Table(names=new_colnames,rows=sub_row)
        table_initialized=True

wd_output_table.add_column('',name='name')
wd_output_table.write(output_filename,delimiter=',',overwrite=True)















