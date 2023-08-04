"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-10-18. Spun-off from select_elBadry_eDR3.py 
on 2023-06-28. 

Take in the el-Badry eDR3 binary catalogue and all of its assorted data and output whatever 
subset of that data with the added info (probably M_G especially).

As is the custom now, this whole thing should be written for Python3 if that wasn't obvious.

The spin-off was done with the intention of implementing quality cuts in the process and 
ignoring elBadry's own classifications of the binaries since he used BP-RP and I intend to use 
G-RP. 

The plan will be to first do a filter on chance alignment probabilities to cut this down to size quite quickly and save myself a headache tryign to have this run on a bunch of binaries that don't matter.

"""


from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.table import join as ATjoin
from astropy.table import hstack as AThstack
from astropy import units as u
from astropy import constants as const
import sys



elBadry_file='all_columns_catalog.fits'

elBadry_dir='/Users/BenKaiser/Desktop/elBadry_catalogue/'

#output_filename='dimWDMS_eDR3.fits'
#output_filename='dimWDMS_eDR3_radec.csv'
#output_filename='dimWDMS_F0toK7_eDR3_highconf.fits'
output_filename='dimWDMS_allMS_minsep4p08_eDR3_highconf.fits'

output_filename=elBadry_dir+output_filename

wd_abs_g_cut= 14.8
#wd_g_rp_cut=1.5#old one as of 2023-06-28
wd_g_rp_cut=1.3

primary_max_bright=10. #brightest that the primary is allowed to be in apparent magnitude.

min_separation=4.08 #in arcseconds
min_separation=min_separation/3600. #converted to degrees as is used in the pairdistance column of the table.

ms_line_points=[
    [1.55,16.75],
    [-0.03,4.58]
    ]

ms_F0K7_abs_g_cut=[3.0,7.4] #inner bounds of the 3 test F0 and K7. Assuming the most extreme of each set of 3 is adequate to keep out the objects outside the range... a bit presumptuous, but such is life.
chance_align_cut=0.1 #the limit set by el-Badry for "high-confidence" binaries. So R<0.1 is high-confidence
wd_val=np.int_(1)
ms_val=np.int_(2) #integer values to be able to indicate what the overall binary is comprised of when summed. I.e. WD+MS binary is 3, WD+WD is 2, MS+MS is 4.
unclassified_val=-5 #default value to fill the column with to ensure any binary with an unclassified component does not mistakenly make it through
wdms_val=wd_val+ms_val #Technically this will really be evaluating if the star is or isn't a white dwarf I suppose because it will be whether the white dwarf is above or below the line that cuts along the H-R diagram
elBadry_file=elBadry_dir+elBadry_file


elBadry_full_table=Table.read(elBadry_file)

highconf=np.where(elBadry_full_table['R_chance_align']<chance_align_cut)
wdms_table=elBadry_full_table[highconf]
print('\n\n\nR<',chance_align_cut)
wdms_table.pprint()

#wdms_indices=np.where(elBadry_full_table['binary_type']=='WDMS') #don't want to use his classifiers


#wdms_table=elBadry_full_table[wdms_indices] #not using his classifiers
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

wdms_table.pprint()

#priming the class integer identifiers

class1_array=np.int_(np.ones(len(wdms_table))*unclassified_val)
class2_array=np.int_(np.ones(len(wdms_table))*unclassified_val)

print('class1_array')
print(class1_array)

wdms_table.add_column(class1_array,name='class1')
wdms_table.add_column(class2_array,name='class2')

def do_quality_cuts(input_table,num=1):
    good_comp=np.where((input_table['parallax_over_error'+str(num)]>10) & (input_table['phot_g_mean_flux_over_error'+str(num)]>10) & (input_table['phot_rp_mean_flux_over_error'+str(num)]>10) & (input_table['g_rp'+str(num)]<1e15) & (input_table['duplicated_source'+str(num)]==False))
    good_table_part=input_table[good_comp]
    
    
    return good_table_part

def classify_targets(input_table, num=1):
    slope=(ms_line_points[1][1]-ms_line_points[0][1])/(ms_line_points[1][0]-ms_line_points[0][0])
    calc_val=slope*(input_table['g_rp'+str(num)]-ms_line_points[0][0])+ms_line_points[0][1]
    print('calc_val',calc_val)
    print(input_table['abs_g_mag'+str(num)])
    wd_indices=np.where(calc_val < input_table['abs_g_mag'+str(num)])
    print('wd_indices',wd_indices)
    input_table['class'+str(num)][wd_indices]=wd_val
    ms_indices=np.where(calc_val >= input_table['abs_g_mag'+str(num)])
    input_table['class'+str(num)][ms_indices]=ms_val
    print('table inside function')
    input_table.pprint()
    return input_table

good_comp1_table=do_quality_cuts(wdms_table,num=1)
print('good_comp1_table')
good_comp1_table.pprint()
print('good_bothtable')
good_bothtable=do_quality_cuts(good_comp1_table,num=2)
good_bothtable.pprint()

#want a minimum separation distance
good_dists=np.where(good_bothtable['pairdistance']>min_separation)
good_bothtable=good_bothtable[good_dists]

print('new table with min separations>',min_separation*3600.,'"')
good_bothtable.pprint()

#sys.exit()


plt.scatter(good_bothtable['g_rp1'],good_bothtable['abs_g_mag1'],s=4)
plt.scatter(good_bothtable['g_rp2'],good_bothtable['abs_g_mag2'],s=4)
plt.plot([ms_line_points[0][0],ms_line_points[1][0]], [ms_line_points[0][1],ms_line_points[1][1]])
plt.ylim(20,0)
plt.show()

plt.scatter(good_bothtable['pairdistance']*3600.,good_bothtable['phot_g_mean_mag2']-good_bothtable['phot_g_mean_mag1'],s=4)
plt.xlabel('separation in arcseconds')
plt.ylabel('gmag2 - gmag1')
plt.show()

#No longer Going to assume the secondary is the white dwarf for our binaries of interest
#assigning integer values to both of the 
good_bothtable=classify_targets(good_bothtable,num=1)
good_bothtable=classify_targets(good_bothtable,num=2)

print('table after the functions run')
good_bothtable.pprint()

#dim_wds=np.where(good_bothtable['abs_g_mag2']>wd_abs_g_cut)
#cool_wdms_table=good_bothtable[dim_wds]

#cool_wdms_table.pprint()

#Now we need to find the overall classification of the binaries
whole_class=good_bothtable['class1']+good_bothtable['class2']


plt.hist(whole_class, bins=np.arange(-6,5,0.5))
plt.show()


wdms_indices=np.where((whole_class>(np.ones(len(good_bothtable))*(wdms_val-0.5)))&(whole_class<(np.ones(len(good_bothtable))*(wdms_val+0.5))))
new_wdms_table=good_bothtable[wdms_indices]

new_wdms_table.pprint()
print(new_wdms_table['binary_type'])

wdsec_ind=np.where(new_wdms_table['class2']==1)
wdsec=new_wdms_table[wdsec_ind]

plt.scatter(new_wdms_table['g_rp1'],new_wdms_table['abs_g_mag1'],s=4)
plt.scatter(new_wdms_table['g_rp2'],new_wdms_table['abs_g_mag2'],s=4)
plt.plot([ms_line_points[0][0],ms_line_points[1][0]], [ms_line_points[0][1],ms_line_points[1][1]])
plt.ylim(20,0)
plt.show()



plt.scatter(wdsec['g_rp1'],wdsec['abs_g_mag1'],s=4)
plt.scatter(wdsec['g_rp2'],wdsec['abs_g_mag2'],s=4)
plt.plot([ms_line_points[0][0],ms_line_points[1][0]], [ms_line_points[0][1],ms_line_points[1][1]])
plt.ylim(20,0)
plt.show()

#sys.exit()

##### Ok now I need to store information indicating which component I believe is the white dwarf. Then in a minute I'll need to do the blue cut below to that set because I seem to have picked up the very tail end of the late main sequence in my "white dwarf" classifications.

def id_wd_component(input_table):
    wd_comp_array=np.int_(np.ones(len(input_table))*unclassified_val)
    primary_wds=np.where(input_table['class1']==wd_val)
    wd_comp_array[primary_wds]=1
    secondary_wds=np.where(input_table['class2']==wd_val)
    wd_comp_array[secondary_wds]=2
    print('wd_comp_array',wd_comp_array)
    input_table.add_column(wd_comp_array,name='wd_comp')
    return input_table

new_wdms_table=id_wd_component(new_wdms_table)
print(new_wdms_table['wd_comp'])


#limit the table to those with white dwarf components that are in the relevant parameter range.
new_wdms_table.add_column(np.int_(np.zeros(len(new_wdms_table))),name='dimwd')
for row in new_wdms_table:
    if ((row['g_rp'+str(row['wd_comp'])]<wd_g_rp_cut) and (row['abs_g_mag'+str(row['wd_comp'])]>wd_abs_g_cut)):
        #print('dim white dwarf')
        row['dimwd']=1
    else:
        pass

dim_wds=np.where(new_wdms_table['dimwd']==1)
new_wdms_table=new_wdms_table[dim_wds]
#print('dim_wds',dim_wds)
print('\n\nWD+MS Binaries with a dim white dwarf as one component\n')
new_wdms_table.pprint()

plt.scatter(new_wdms_table['g_rp1'],new_wdms_table['abs_g_mag1'],s=4)
plt.scatter(new_wdms_table['g_rp2'],new_wdms_table['abs_g_mag2'],s=4)
plt.plot([ms_line_points[0][0],ms_line_points[1][0]], [ms_line_points[0][1],ms_line_points[1][1]])
plt.ylim(20,0)
plt.show()


##############################
"""
Ok so we have now limited our sample to those binaries that are likely to actually host a dim white dwarf in addition to a main sequence (or potentially actually red giant) component.
The next step would be to come up with some function that will eliminate binary pairs that have overly bright main sequence companions for the separation distance of the two stars.
I imagine the brightness falloff is something exponential most likely given that it's roughly the edge of a Gaussian (and then ya know the whole Airy disk Legendre's), but 
Up to now there are 388 results that get returned.

I'll probably need to finally do a requery in Gaia to check for crowded fields and unrelated nearby contaminants that could be screwing up colors.
I suspect I'm going to lose a ton of these objects in the coming days.
"""


###############################
sys.exit()



#blue_wds=np.where(cool_wdms_table['g_rp2']<wd_g_rp_cut)
#cool_wdms_table=cool_wdms_table[blue_wds]

#should remove the objects whose primary is actually a WD.

#ms_primaries=np.where(cool_wdms_table['abs_g_mag1']<((17.77-7.5)/(1.3))*cool_wdms_table['g_rp1']+7.5)

#cool_wdms_table=cool_wdms_table[ms_primaries]
#print('MS primaries only:')
#cool_wdms_table.pprint()


#plt.scatter(cool_wdms_table['g_rp1'],cool_wdms_table['abs_g_mag1'],s=4)
#plt.scatter(cool_wdms_table['g_rp2'],cool_wdms_table['abs_g_mag2'],s=4)
#plt.ylim(20,0)
#plt.show()

#preK7_primaries=np.where((cool_wdms_table['abs_g_mag1']< ms_F0K7_abs_g_cut[1]))
#cool_wdms_table=cool_wdms_table[preK7_primaries]

#postF0_primaries=np.where((cool_wdms_table['abs_g_mag1']> ms_F0K7_abs_g_cut[0]))
#cool_wdms_table=cool_wdms_table[postF0_primaries]
#print('\n\n\nF0 to K7 primaries only')
#cool_wdms_table.pprint()


#plt.scatter(cool_wdms_table['g_rp1'],cool_wdms_table['abs_g_mag1'],s=4)
#plt.scatter(cool_wdms_table['g_rp2'],cool_wdms_table['abs_g_mag2'],s=4)
#plt.ylim(20,0)
#plt.show()

#highconf=np.where(cool_wdms_table['R_chance_align']<chance_align_cut)
#cool_wdms_table=cool_wdms_table[highconf]
#print('\n\n\nR<',chance_align_cut)
#cool_wdms_table.pprint()

sys.exit()

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
