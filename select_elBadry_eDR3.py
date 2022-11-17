"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-10-18

Take in the el-Badry eDR3 binary catalogue and all of its assorted data and output whatever 
subset of that data with the added info (probably M_G especially).

As is the custom now, this whole thing should be written for Python3 if that wasn't obvious.


"""


from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.table import join as ATjoin
from astropy.table import hstack as AThstack
from astropy import units as u
from astropy import constants as const



elBadry_file='all_columns_catalog.fits'

elBadry_dir='/Users/BenKaiser/Desktop/elBadry_catalogue/'

#output_filename='dimWDMS_eDR3.fits'
#output_filename='dimWDMS_eDR3_radec.csv'
output_filename='dimWDMS_F0toK7_eDR3_highconf.fits'

output_filename=elBadry_dir+output_filename

wd_abs_g_cut= 14.8
wd_g_rp_cut=1.5

ms_F0K7_abs_g_cut=[3.0,7.4] #inner bounds of the 3 test F0 and K7. Assuming the most extreme of each set of 3 is adequate to keep out the objects outside the range... a bit presumptuous, but such is life.
chance_align_cut=0.1 #the limit set by el-Badry for "high-confidence" binaries. So R<0.1 is high-confidence

elBadry_file=elBadry_dir+elBadry_file


elBadry_full_table=Table.read(elBadry_file)


wdms_indices=np.where(elBadry_full_table['binary_type']=='WDMS')


wdms_table=elBadry_full_table[wdms_indices]
#wdms_table.pprint()
print(len(wdms_table))

def get_abs_mag(band='g',comp=1):
    #abs_mag=wdms_table['phot_'+band+'_mean_mag'+str(comp)]+5*np.log10(wdms_table['parallax'+str(comp)])-10
    abs_mag=wdms_table['phot_'+band+'_mean_mag'+str(comp)]+5*np.log10(wdms_table['parallax'+'1'])-10 #El Badry uses the parallax of the primary to get the absolute magnitude of the secondary because it usually has better precision I'm pretty sure...
    return abs_mag

abs_g_mag1=get_abs_mag(comp=1)
abs_g_mag2=get_abs_mag(comp=2)

wdms_table.add_column(abs_g_mag1,name='abs_g_mag1')
wdms_table.add_column(abs_g_mag2,name='abs_g_mag2')


#Going to assume the secondary is the white dwarf for our binaries of interest

dim_wds=np.where(abs_g_mag2>wd_abs_g_cut)
cool_wdms_table=wdms_table[dim_wds]

#cool_wdms_table.pprint()


blue_wds=np.where(cool_wdms_table['g_rp2']<wd_g_rp_cut)
cool_wdms_table=cool_wdms_table[blue_wds]

#should remove the objects whose primary is actually a WD.

ms_primaries=np.where(cool_wdms_table['abs_g_mag1']<((17.77-7.5)/(1.3))*cool_wdms_table['g_rp1']+7.5)

cool_wdms_table=cool_wdms_table[ms_primaries]
print('MS primaries only:')
cool_wdms_table.pprint()

preK7_primaries=np.where((cool_wdms_table['abs_g_mag1']< ms_F0K7_abs_g_cut[1]))
cool_wdms_table=cool_wdms_table[preK7_primaries]

postF0_primaries=np.where((cool_wdms_table['abs_g_mag1']> ms_F0K7_abs_g_cut[0]))
cool_wdms_table=cool_wdms_table[postF0_primaries]
print('\n\n\nF0 to K7 primaries only')
cool_wdms_table.pprint()

highconf=np.where(cool_wdms_table['R_chance_align']<chance_align_cut)
cool_wdms_table=cool_wdms_table[highconf]
print('\n\n\nR<',chance_align_cut)
cool_wdms_table.pprint()

cool_wdms_table.write(output_filename)

secondary_table=Table([cool_wdms_table['ra2'],cool_wdms_table['dec2']],names=['ra','dec'])

secondary_table.write(output_filename)
#cool_wdms_table.header()

new_table_list=[]
for name in cool_wdms_table.col():
    print(name)
#cool_wdms_table.write('el_badry_dimWDMS.csv')


plt.scatter(cool_wdms_table['g_rp1'],cool_wdms_table['abs_g_mag1'],s=4)
plt.scatter(cool_wdms_table['g_rp2'],cool_wdms_table['abs_g_mag2'],s=4)
plt.ylim(20,0)
plt.show()


plt.scatter(cool_wdms_table['g_rp2'],cool_wdms_table['phot_g_mean_mag2'],s=4)
plt.show()

#plt.scatter(abs_g_mag1,abs_g_mag2)
#plt.xlabel('abs_g_mag1')
#plt.ylabel('abs_g_mag2')
#plt.show()

#plt.scatter(abs_g_mag1,abs_g_mag2-abs_g_mag1,s=4)
#plt.xlabel('abs_g_mag1')
##plt.ylabel('abs_g_mag2')
#plt.show()
