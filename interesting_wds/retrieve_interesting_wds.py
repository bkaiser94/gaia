"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-07-30

This file should load the table I made with retrieve_faint_wds.py and then plot the points and filter by the cuts from the MORDOR survey riy cuts.



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
from astropy.table import Table, Column, vstack, join, MaskedColumn
#import scipy.interpolate as scinterp
import time
#import pyvo
#from astroquery.vizier import Vizier

#from astroquery.gaia import Gaia
#from astroquery.xmatch import XMatch


sys.path.append('../')
sys.path.append('/Users/BenKaiser/Desktop/radial_velocity_calculations/')
import cal_params as cp
#import spec_plot_tools as spt

plt.rc('lines',linewidth=0.5)

input_table=Table.read('gf21_GaiaeDR3_faintWDs_gaiaadded.csv')


def wddm_cut(ri):
    return 0.03*ri+0.52


def interesting_cut(ri):
    return 1.66*ri-0.53

def interesting_50(ri):
    return 1.66*ri-0.93

def get_mag(table,band,survey='ps1'):
    if survey=='ps1':
        mags=table['ps1_'+band+'_mean_psf_mag']
        errs=table['ps1_'+band+'_mean_psf_mag_error']
    elif survey=='sdss':
        mags=table[band+'mag']
        errs=table[band+'mag']*0

    return mags, errs

def get_color(table,colorlist):
        mag_list=[]
        err_list=[]
        for item in colorlist:
            #print('item',item,item[0])
            mag, err=get_mag(table,item[0],survey=item[1])
            mag_list.append(mag)
            err_list.append(err)
        np.array(mag_list), np.array(err_list)
        color=mag_list[0]-mag_list[1]
        color_err=np.sqrt(err_list[0]**2+err_list[1]**2)
        return color, color_err

def plot_color_color(table,colorlist1,colorlist2,plot_color='k',label=''):
    """
    From gaia/mordor_survey_code/grizy_spectype_cuts.ipynb
    plot band1-band2 on x , band3-band4 on y
    colorlist1: [[band1,survey1],[band2,survey2]]
    colorlist2:[[band3,survey3],[band4,survey4]]
    """
    def get_color(colorlist):
        mag_list=[]
        err_list=[]
        for item in colorlist:
            #print('item',item,item[0])
            mag, err=get_mag(table,item[0],survey=item[1])
            mag_list.append(mag)
            err_list.append(err)
        return mag_list, err_list
    mag_list1, err_list1=get_color(colorlist1)
    mag_list2, err_list2=get_color(colorlist2)
    #print('err_list1',err_list1)
    #print('err_list2',err_list2)
    max_err=np.argmax(np.sqrt(err_list1[0]**2+err_list1[1]**2))
    max_err2=np.argmax(np.sqrt(err_list2[0]**2+err_list2[1]**2))
    print(((mag_list2[0]-mag_list2[1])+3*np.sqrt(err_list2[0]**2+err_list2[1]**2))[max_err2])

    plt.errorbar(mag_list1[0]-mag_list1[1],mag_list2[0]-mag_list2[1],xerr= np.sqrt(err_list1[0]**2+err_list1[1]**2), yerr=np.sqrt(err_list2[0]**2+err_list2[1]**2) ,color=plot_color,label=label, linestyle='None',marker='o',markersize=2,alpha=0.3,markeredgecolor='none')
    
    #plt.errorbar(mag_list1[0]-mag_list1[1],mag_list2[0]-mag_list2[1] ,color=plot_color,label=label, linestyle='None',marker='o',markersize=2,alpha=0.1,markeredgecolor='none')
    
    def make_axis_label(colorlist):
        axis_label=colorlist[0][1]+'_'+colorlist[0][0]+' - '+colorlist[1][1]+'_'+colorlist[1][0]
        return axis_label
    #for x,y,s in zip(mag_list1[0]-mag_list1[1],mag_list2[0]-mag_list2[1],table['name']):
        #plt.text(x,y,s,color=plot_color)
    #plt.text(mag_list1[0]-mag_list1[1],mag_list2[0]-mag_list2[1],table['name'],color=plot_color)
    stringx=make_axis_label(colorlist1)
    stringy=make_axis_label(colorlist2)
    plt.xlabel(stringx)
    plt.ylabel(stringy)

    
    return

list1=[['r','ps1'],['i','ps1']]
list2=[['i','ps1'],['y','ps1']]

ri,ri_err=get_color(input_table,list1)
iy,iy_err=get_color(input_table,list2)


print('type(ri)',type(ri))
print('ri',ri)
ri_col=MaskedColumn(ri,name='ps1_r_i')
print('ri_col',ri_col)


iy_col=MaskedColumn(iy,name='ps1_i_y')
iy_err_col=MaskedColumn(iy_err,name='ps1_i_y_err')
ri_err_col=MaskedColumn(ri_err,name='ps1_r_i_err')
print('len(ri_col)',len(ri_col))
print('len(input_table)',len(input_table))



input_table.add_columns([ri_col,ri_err_col,iy_col,iy_err_col])
ri_err_no_mask=~input_table['ps1_r_i_err'].mask
iy_err_no_mask=~input_table['ps1_i_y_err'].mask

print('ri_err_no_mask',ri_err_no_mask,np.sum(ri_err_no_mask))

no_mask=ri_err_no_mask*iy_err_no_mask

good_table=input_table[no_mask]

interesting_wd_indices=np.where((good_table['ps1_i_y']<wddm_cut(good_table['ps1_r_i'])) & (good_table['ps1_i_y']<interesting_cut(good_table['ps1_r_i']))& (good_table['ps1_i_y']>=interesting_50(good_table['ps1_r_i'])))
interesting_table=good_table[interesting_wd_indices]

print('interesting table')
interesting_table.pprint()
print(len(interesting_table))

vinterest_wd_indices=np.where((good_table['ps1_i_y']<wddm_cut(good_table['ps1_r_i'])) & (good_table['ps1_i_y']<interesting_50(good_table['ps1_r_i'])))
vinterest_table=good_table[vinterest_wd_indices]
print('Very interesting table')
vinterest_table.pprint()
print(len(vinterest_table))
xvals=np.array([-4.5,4.5])

print('start min maxes')
print(np.nanmin(vinterest_table['ps1_r_i']))
print(np.nanmax(vinterest_table['ps1_r_i']))
print(np.nanmin(vinterest_table['ps1_i_y']))
print(np.nanmax(vinterest_table['ps1_i_y']))
print('end of min maxes')
#plot_color_color(input_table,list1,list2)
plot_color_color(good_table,list1,list2)
plot_color_color(interesting_table,list1,list2,label=str(len(interesting_table))+' pretty interesting wds',plot_color='blue')
plot_color_color(vinterest_table,list1,list2,label=str(len(vinterest_table))+' very interesting wds',plot_color='magenta')
wddm_yvals=wddm_cut(xvals)
interesting_yvals=interesting_cut(xvals)
interesting_50_yvals=interesting_50(xvals)
plt.plot(xvals,wddm_yvals,color=cp.sp_type_colors['WDdM'],label='WD+dM cut')
plt.plot(xvals,interesting_yvals, color=cp.sp_type_colors['DC'],label='Interesting WD cut')
plt.plot(xvals,interesting_50_yvals, color=cp.sp_type_colors['DC'],linestyle='--')
plt.legend()
#plt.xlim(-4.5,4.5)
#plt.ylim(-3,5.3)
plt.xlim(-1,3)
plt.ylim(-1.5,2.5)
plt.show()



#vinterest_table.write('gf21_GaiaeDR3_faintWDs_gaiaadded_veryinterestingwds.csv')
#interesting_table.write('gf21_GaiaeDR3_faintWDs_gaiaadded_prettyinterestingwds.csv')




