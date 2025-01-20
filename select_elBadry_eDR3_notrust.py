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
from astropy.table import vstack as ATvstack
from astropy import units as u
from astropy import constants as const
import sys



elBadry_file='all_columns_catalog.fits'

elBadry_dir='/Users/BenKaiser/Desktop/elBadry_catalogue/'

#output_filename='dimWDMS_eDR3.fits'
#output_filename='dimWDMS_eDR3_radec.csv'
#output_filename='dimWDMS_F0toK7_eDR3_highconf.fits'
#output_filename='dimWDMS_allMS_minsep4p08_eDR3_highconf.fits'
#output_filename='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust.fits'
#output_filename='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWDradec.csv'
#output_filename='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWDradec.csv'
#output_filename='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD.fits'
output_filename='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_loosecuts.fits'


output_filename=elBadry_dir+output_filename

wd_abs_g_cut= 14.8
#wd_g_rp_cut=1.5#old one as of 2023-06-28
wd_g_rp_cut=1.3

#primary_max_bright=10. #brightest that the primary is allowed to be in apparent magnitude.

#min_separation=4.08 #in arcseconds
#min_separation=min_separation/3600. #converted to degrees as is used in the pairdistance column of the table.

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


#################### binary separation function values#########

close_mag_diff=-4.5 #constant value inside the closest radius (these are basically completely incapable of being disentangled in the reduction (or rather removing the contaminant flux so the contaminant flux needs to be much less than the target. Hence the negative value
close_mag_radius=4. #in arcseconds
mid_mag_diff=4.6 
mid_mag_out_radius=13. #in arcseconds

mag_coeff=-9.8
mag_center=7.
mag_spread=15.
mag_offset=14.6




#################################################


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
    #good_comp=np.where((input_table['parallax_over_error'+str(num)]>10) & (input_table['phot_g_mean_flux_over_error'+str(num)]>10) & (input_table['phot_rp_mean_flux_over_error'+str(num)]>10) & (input_table['g_rp'+str(num)]<1e15) & (input_table['duplicated_source'+str(num)]==False)) #original cuts
    good_comp=np.where((input_table['parallax_over_error'+str(num)]>5) & (input_table['phot_g_mean_flux_over_error'+str(num)]>5) & (input_table['phot_rp_mean_flux_over_error'+str(num)]>5) & (input_table['g_rp'+str(num)]<1e15) & (input_table['duplicated_source'+str(num)]==False)) #loose cuts to be able to recover analogues to WDJ0356-2255 and SDSS J1330+6435 with their broad Na I D features and poorer signal to noise. Implemented 2023-10-05
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
#I'm removing this because it's not a minimum separation distance we need (in a sense) it is a minimum magnitude difference for various separations
#good_dists=np.where(good_bothtable['pairdistance']>min_separation)
#good_bothtable=good_bothtable[good_dists]

#print('new table with min separations>',min_separation*3600.,'"')
#good_bothtable.pprint()

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
new_wdms_table.add_column(np.float_(np.zeros(len(new_wdms_table))),name='mag_diff')
for row in new_wdms_table:
    if row['wd_comp']==1:
        other_comp=2
    else:
        other_comp=1
    mag_diff=row['phot_g_mean_mag'+str(row['wd_comp'])]-row['phot_g_mean_mag'+str(other_comp)]
    row['mag_diff']=mag_diff
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
distance_arcseconds=new_wdms_table['pairdistance']*3600.
close_pairs=np.where(distance_arcseconds<close_mag_radius)
close_good_inds=np.where(new_wdms_table[close_pairs]['mag_diff']<close_mag_diff)
close_good_pairs=new_wdms_table[close_pairs][close_good_inds]

mid_pairs=np.where((distance_arcseconds>= close_mag_radius) & (distance_arcseconds<=mid_mag_out_radius))
mid_good_inds=np.where(new_wdms_table[mid_pairs]['mag_diff']<mid_mag_diff)
mid_good_pairs=new_wdms_table[mid_pairs][mid_good_inds]

def far_mag_func(dist_column):
    mag_boundary=mag_coeff*np.exp(-1.*(dist_column-mag_center)**2./(2*mag_spread**2.))+mag_offset
    return mag_boundary

far_pairs=np.where(distance_arcseconds>mid_mag_out_radius)
far_good_inds=np.where(new_wdms_table[far_pairs]['mag_diff']<far_mag_func(distance_arcseconds[far_pairs]))
far_good_pairs=new_wdms_table[far_pairs][far_good_inds]

print('close_good_pairs')
close_good_pairs.pprint()

print('\n\nmid_good_pairs')
mid_good_pairs.pprint()

print('\n\nfar_good_pairs')
far_good_pairs.pprint()

markersize=6
plt.scatter(new_wdms_table['pairdistance']*3600., new_wdms_table['mag_diff'], label='full WDMS table before magdiff cuts',s=markersize)
plt.scatter(close_good_pairs['pairdistance']*3600., close_good_pairs['mag_diff'],label='close pairs',s=markersize)
plt.scatter(mid_good_pairs['pairdistance']*3600., mid_good_pairs['mag_diff'],label='mid pairs',s=markersize)
plt.scatter(far_good_pairs['pairdistance']*3600., far_good_pairs['mag_diff'], label='far pairs',s=markersize)

close_test=np.linspace(0,close_mag_radius,100)
plt.plot(close_test, np.ones(close_test.shape)*close_mag_diff)
mid_test=np.linspace(close_mag_radius,mid_mag_out_radius,100)
mid_vals=np.ones(mid_test.shape)*mid_mag_diff
plt.plot(mid_test,mid_vals)
far_test=np.linspace(mid_mag_out_radius,400.,1000)
far_vals=far_mag_func(far_test)
plt.plot(far_test,far_vals)

plt.plot(31.,9.41,marker='*',color='k',label='SDSS J1312-0229 and its companion',markersize=10)

plt.xlabel('Separation (arcseconds)')
plt.ylabel('G Mag difference (WD - Companion)')
plt.title('Attempt at Magnitude Difference Cut based on Separation of the Binary Components')
plt.legend()
plt.show()

plt.scatter(new_wdms_table['l1'],new_wdms_table['b1'],label='new_wdms_table',s=markersize)
plt.scatter(close_good_pairs['l1'],close_good_pairs['b1'],label='close_good_pairs',s=markersize)
plt.scatter(mid_good_pairs['l1'],mid_good_pairs['b1'],label='mid_good_pairs',s=markersize)
plt.scatter(far_good_pairs['l1'],far_good_pairs['b1'],label='far_good_pairs',s=markersize)

plt.title('Galactic Coords of Selected Targets')
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.legend()
plt.show()

###############################
#sys.exit()

output_fulltable=ATvstack([close_good_pairs,mid_good_pairs,far_good_pairs])

"""
Now we need to make a table to output to use for future efforts. This table should presumably only include the
white dwarf component... but maybe it would be helpful to retain the full table for plotting purposes. I suppose I could have two different outputs
but I only run it with one active at a time to prevent shenanigans and mistakes.

So first I guess I'll get the version that doesn't trim down to just the white dwarfs running.

"""

print('\n\noutput_fulltable')
output_fulltable.pprint()
#output_fulltable.write(output_filename,overwrite=True)

sys.exit()
"""
Ok now I need to get the version that only contains the white dwarfs running... 

so I'll probably have to do this row-by-row and have it check the wd_comp column each time for what values it should retain. I guess I need to make 
a table in advance to have this populate though... But that could just be a copy of the original table. But that actually won't work because that is too 
many columns. I need to make a version of the table to populate that doesn't have double the columns (while still retaining the joint binary 
information maybe...)

I actually only want to have this output the WD RA and DEC because the point of only retaining the white dwarfs is to be able to construct the target 
list and more importantly requery into Gaia DR3 to make sure it is still consistent with the minimum separation and magnitude function we used for 
the binary companions. Also probably going to use the stellar density cut still because I doubt other observers will be able to target efficiently in 
crowded fields. Oh limitations...

"""
#full_colnames=output_fulltable.colnames()
#print('colnames:',full_colnames)
#new_colnames=[]
#for col in full_colnames:
    #if col[-1]=='1':
        #new_colnames.append(col[:-1])
        ##this method works because there is no DR1 indicator. Had I used the secondary component's columns this would have accidentally picked up all of the DR2 columns. I also had it exclude the last character because that is usually the number itself.
    #else:
        #pass



#output_copy=output_fulltable.copy()
#print('\nTable copied\n')

sys.exit()

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

##########################################


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
