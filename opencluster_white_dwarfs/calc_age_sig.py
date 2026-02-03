"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-09-26 from the plotting code

Calculate the significance in multiples of sigma of the difference in the cooling age of the white dwarf and the cluster age.

Then output that into a new file as an additional column of the table.



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
from astropy.table import Table, Column, vstack, join, hstack
#import scipy.interpolate as scinterp
import time
import pyvo
import cluster_utilities as cu


sys.path.append('../')


#input_file='HR24members_GF21maincat_crossmatch_simbadadded_wdagesadded.csv'
input_file='HR24members_GF21maincat_crossmatch_simbadadded_mwddadded_wdagesadded.csv' #MWDD parameter-based ages,
#input_file='HR24members_GF21maincat_crossmatch_simbadadded_mwddadded_wdagesaddedonlyDAmodels.csv'

input_table=Table.read(input_file)

#input_table=input_table[:5]

input_table['total_age_median'].pprint()
print(input_table['total_age_median']-(input_table['total_age_err_low']))
print(np.log10(input_table['total_age_median']),np.log10(input_table['total_age_median']-(input_table['total_age_err_low'])))
print(np.log10(input_table['total_age_median'])+9)
print(np.log10(input_table['total_age_median']/(input_table['total_age_median']-input_table['total_age_err_low']))+9)
print(np.log10(input_table['total_age_median'])-np.log10(input_table['total_age_median']-input_table['total_age_err_low'])+9)



def err_quad_sum(table):
    
    """
    Assuming the white dwarf cooling age is greater than the cluster age because that is the scenario we care about, so we're going to add the upper cluster age error bar to the lower cooling age error bar in quadrature
    
    """
    linear_cluster_up_err=10**table['logAge84']-10**table['logAge50'] #cluster upper error in years
    cooling_err=table['cooling_age_err_low']*1e9 #converting lower cooling age error to years
    
    return np.sqrt(linear_cluster_up_err**2+cooling_err**2)



def calc_age_sig(table):
    quad_err=err_quad_sum(table) #error in years on the difference in the cooling age and cluster age... assuming the cooling age is greater than the cluster age technically... it's probably OK for doing lesser cooling ages, but technically not right... assuming I'm right for the other scenario, which is admittedly not a guarantee.
    
    difference=table['cooling_age_median']*1E9-10**table['logAge50'] #difference between the cooling age and the cluster age in years
    
    
    
    
    return difference*1E-9, quad_err*1E-9, difference/quad_err



if __name__ == '__main__':
    
    age_diff, age_diff_err,age_diff_sig=calc_age_sig(input_table)
    input_table['age_diff']=age_diff
    input_table['age_diff_err']=age_diff_err
    input_table['age_diff_sig']=age_diff_sig
    output_name=input_file.split('.')[0]+'_agediffsig.csv'
    input_table.write(output_name)





