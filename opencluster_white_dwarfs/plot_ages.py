"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-09-17

Plot the cluster ages and white dwarf total/cooling ages together to hopefully reveal outlier points.



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
input_file='HR24members_GF21maincat_crossmatch_simbadadded_mwddadded_wdagesadded_agediffsig.csv' #MWDD parameter-based ages,
#input_file='HR24members_GF21maincat_crossmatch_simbadadded_mwddadded_wdagesaddedonlyDAmodels.csv'

prob_threshold=0.999


input_table=Table.read(input_file)

#input_table=input_table[:5]

input_table['total_age_median'].pprint()
print(input_table['total_age_median']-(input_table['total_age_err_low']))
print(np.log10(input_table['total_age_median']),np.log10(input_table['total_age_median']-(input_table['total_age_err_low'])))
print(np.log10(input_table['total_age_median'])+9)
print(np.log10(input_table['total_age_median']/(input_table['total_age_median']-input_table['total_age_err_low']))+9)
print(np.log10(input_table['total_age_median'])-np.log10(input_table['total_age_median']-input_table['total_age_err_low'])+9)


#sig_havers_inds=np.where(input_table['age_diff_sig']>-100)
sig_havers_inds=input_table['age_diff_sig'].mask
sig_havers=input_table[~sig_havers_inds]
print('\n\nsig_havers')
sig_havers['age_diff_sig'].pprint()

high_prob_inds=np.where(sig_havers['Prob']>=prob_threshold)
high_prob=sig_havers[high_prob_inds]

impossible_inds=np.where(high_prob['age_diff_sig']>=3)
impossible_table=high_prob[impossible_inds]

print('Num impossibly cool WDs with >='+str(prob_threshold)+ ' membership prob',len(impossible_table))
print('Num of WDs with >=' + str(prob_threshold)+' membership and age_diff_sig values',len(high_prob))
print('Fraction of >=' +str(prob_threshold)+ ' membership that have incompatible *cooling ages* (so even more signficant with whatever main-sequence contribution',len(impossible_table)/len(high_prob))


def plot_ages(input_table=input_table,age_string='total',marker='o',label='default',alpha=0.4):
    
    #plt.errorbar(input_table['logAge50'],input_table['total_age_median'],xerr=[input_table['logAge50']-input_table['logAge16'],input_table['logAge84']-input_table['logAge50']], yerr=[input_table['total_age_err_low'],input_table['total_age_err_high']],linestyle='None',marker='o')
    if label=='default':
        label=age_string+' Age'
    else:
        label=label
    plt.errorbar(input_table['logAge50'],np.log10(input_table[age_string+'_age_median'])+9,xerr=[input_table['logAge50']-input_table['logAge16'],input_table['logAge84']-input_table['logAge50']], yerr=[np.log10(input_table[age_string+'_age_median']/(input_table[age_string+'_age_median']-input_table[age_string+'_age_err_low'])),np.log10((input_table[age_string+'_age_median']+input_table[age_string+'_age_err_high'])/input_table[age_string+'_age_median'])],linestyle='None',label=label,marker=marker,alpha=alpha)
    return

def plot_age_v_Teff(input_table=input_table,age_string='total',marker='o',label='default',alpha=0.4, log=True):
    
    #plt.errorbar(input_table['logAge50'],input_table['total_age_median'],xerr=[input_table['logAge50']-input_table['logAge16'],input_table['logAge84']-input_table['logAge50']], yerr=[input_table['total_age_err_low'],input_table['total_age_err_high']],linestyle='None',marker='o')
    if label=='default':
        label=age_string+' Age'
    else:
        label=label
    if log:
        plt.errorbar(input_table['mwdd_teff'],np.log10(input_table[age_string+'_age_median'])+9,xerr=input_table['mwdd_Dteff'], yerr=[np.log10(input_table[age_string+'_age_median']/(input_table[age_string+'_age_median']-input_table[age_string+'_age_err_low'])),np.log10((input_table[age_string+'_age_median']+input_table[age_string+'_age_err_high'])/input_table[age_string+'_age_median'])],linestyle='None',label=label,marker=marker,alpha=alpha)
    else:
        plt.errorbar(input_table['mwdd_teff'],input_table[age_string+'_age_median'],xerr=input_table['mwdd_Dteff'], yerr=[input_table[age_string+'_age_err_low'],input_table[age_string+'_age_err_high']],linestyle='None',label=label,marker=marker,alpha=alpha)
    return


def plot_agesig_v_clusterage(input_table=input_table,marker='o',alpha=1, log=True,label=''):
    print('type(input_table)',type(input_table))
    plt.errorbar(input_table['logAge50'],input_table['age_diff_sig'],xerr=[input_table['logAge50']-input_table['logAge16'],input_table['logAge84']-input_table['logAge50']],linestyle='None',label=label,marker=marker,alpha=alpha)
    #for row in input_table:
        #plt.text(row['logAge50'],row['age_diff_sig'],row['Name'])
    plt.xlabel('log10(Cluster Age)')
    plt.ylabel('Cooling Age Discrepancy Significance\n[(cooling age- cluster age)/(cooling_err^2+cluster_err^2)^(1/2)]')
    return

def plot_agesig_v_prob(input_table=input_table,marker='o',alpha=1, log=True,label=''):
    print('type(input_table)',type(input_table))
    plt.errorbar(input_table['Prob'],input_table['age_diff_sig'],linestyle='None',label=label,marker=marker,alpha=alpha)
    for row in input_table:
        plt.text(row['Prob'],row['age_diff_sig'],row['Name'])
    plt.xlabel('Probability of WD being a member of association per Hunt and Reffert (2024)')
    plt.ylabel('Cooling Age Discrepancy Significance\n[(cooling age- cluster age)/(cooling_err^2+cluster_err^2)^(1/2)]')
    return



rsg5=np.where(input_table['Name']=='RSG_5')


dahs=cu.limit_spectype(input_table,spectype='DAH',flexible_spectype=False,spectype_colname='mwdd_spectype')
nonDAs=cu.limit_spectype(input_table,spectype='non-DA',flexible_spectype=False,spectype_colname='wdwarfdate_model')
DAs=cu.limit_spectype(input_table,spectype='DA',flexible_spectype=True,spectype_colname='mwdd_spectype')


plot_agesig_v_clusterage(input_table=DAs,label='DA',alpha=0.4)
plot_agesig_v_clusterage(input_table=nonDAs,label='non-DA',marker='D',alpha=0.4)
plot_agesig_v_clusterage(input_table=input_table[rsg5],label='RSG_5',marker='H')
plot_agesig_v_clusterage(input_table=dahs,label='DAH',marker='s')
plt.axhline(y=3,linestyle='--',color='k')
plt.text(9,3,'3-sigma')
plt.legend()
plt.show()

plot_agesig_v_prob(input_table=DAs,label='DA',alpha=0.4)
plot_agesig_v_prob(input_table=nonDAs,label='non-DA',marker='D',alpha=0.4)
plot_agesig_v_prob(input_table=input_table[rsg5],label='RSG_5',marker='H')
plot_agesig_v_prob(input_table=dahs,label='DAH',marker='s')
plt.axhline(y=3,linestyle='--',color='k')
plt.text(0.55,3,'3-sigma')
plt.legend()
plt.show()




plot_age_v_Teff(input_table=DAs,age_string='cooling',label='All DA white dwarfs')
plot_age_v_Teff(input_table=nonDAs,age_string='cooling',label='All non-DA white dwarfs',marker='D')
plot_age_v_Teff(input_table=input_table[rsg5],age_string='cooling',marker='H',label='RSG-5',alpha=1)
plot_age_v_Teff(input_table=dahs,age_string='cooling',marker='s',label='DAH',alpha=1)
plt.xlabel('Teff (K)')
plt.ylabel('Log10(WD Cooling Age) [yrs]')
plt.title(input_file)
plt.legend()
plt.show()

plot_age_v_Teff(input_table=DAs,age_string='cooling',label='All DA white dwarfs',log=False)
plot_age_v_Teff(input_table=nonDAs,age_string='cooling',label='All non-DA white dwarfs',marker='D',log=False)
plot_age_v_Teff(input_table=input_table[rsg5],age_string='cooling',marker='H',label='RSG-5',alpha=1,log=False)
plot_age_v_Teff(input_table=dahs,age_string='cooling',marker='s',label='DAH',alpha=1,log=False)
plt.xlabel('Teff (K)')
plt.ylabel('WD Cooling Age (Gyrs)')
plt.title(input_file)
plt.legend()
plt.show()




plt.xlabel('Log10(Cluster Age) [yrs]')
plt.ylabel('Log10(WD Cooling Age) [yrs]')
plt.plot([0,10],[0,10],linestyle='--',color='k')
#plot_ages()
plot_ages(input_table=DAs,age_string='cooling',label='All DA white dwarfs')
plot_ages(input_table=nonDAs,age_string='cooling',label='All non-DA white dwarfs',marker='D')
plot_ages(input_table=input_table[rsg5],age_string='cooling',marker='H',label='RSG-5',alpha=1)
plot_ages(input_table=dahs,age_string='cooling',marker='s',label='DAH',alpha=1)
plt.legend()
plt.text(6.1,9.5,'White Dwarfs above the dashed line\nare too old for their open clusters\nbased on cooling age alone.')
plt.xlim(6,10)
plt.ylim(6,10)
plt.title(input_file)
plt.show()

plt.xlabel('Log10(Cluster Age) [yrs]')
plt.ylabel('Log10(WD Total Age) [yrs]')
plt.plot([0,11],[0,11],linestyle='--',color='k')
#plot_ages()
plot_ages(input_table=DAs,age_string='total',label='All DA white dwarfs')
plot_ages(input_table=nonDAs,age_string='total',label='All non-DA white dwarfs',marker='D')
plot_ages(input_table=input_table[rsg5],age_string='total',marker='H',label='RSG-5',alpha=1)
plot_ages(input_table=dahs,age_string='total',marker='s',label='DAH',alpha=1)
plt.legend()
plt.text(6.1,9.5,'White Dwarfs above the dashed line\nare too old for their open clusters\nbased on total age.')
plt.xlim(6,11)
plt.ylim(6,11)
plt.title(input_file)
plt.show()




dahs=np.where(input_table['simbad_sp_type']=='DAH')
plt.xlabel('Log10(Cluster Age) [yrs]')
plt.ylabel('Log10(WD Cooling Age) [yrs]')
plt.plot([0,10],[0,10],linestyle='--',color='k')
#plot_ages()
plot_ages(age_string='cooling',label='All MWDD-parameter white dwarfs')
plot_ages(input_table=input_table[rsg5],age_string='cooling',marker='*',label='RSG-5',alpha=1)
plot_ages(input_table=input_table[dahs],age_string='cooling',marker='s',label='DAH',alpha=1)
plt.legend()
plt.text(6.1,9.5,'White Dwarfs in this triangle\nare too old for their open clusters\nbased on cooling age alone.')
plt.xlim(6,10)
plt.ylim(6,10)
plt.title(input_file)
plt.show()


plt.xlabel('Log10(Cluster Age) [yrs]')
plt.ylabel('Log10(WD Total Age) [yrs]')
plt.plot([0,10],[0,10],linestyle='--',color='k')
#plot_ages()
plot_ages(age_string='cooling',label='All MWDD-parameter white dwarfs')
plot_ages(input_table=input_table[rsg5],age_string='cooling',marker='*',label='RSG-5',alpha=1)
plot_ages(input_table=input_table[dahs],age_string='cooling',marker='s',label='DAH',alpha=1)
plt.legend()
plt.text(6.1,9.5,'White Dwarfs in this triangle\nare too old for their open clusters\nbased on Total age.')
plt.xlim(6,10)
plt.ylim(6,10)
plt.title(input_file)
plt.show()
