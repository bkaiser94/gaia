"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-09-17


Plot a given Gaia target list as a scatter plot with error bars for a variety of parameters. Will also eventually include 
quality cuts overlaid perhaps, and maybe a background plotting version.

Basically a hybrid of plot_alt_cmd.py and plot_data_checks.py

Actually this might end up being like plot_alt_cmd.py essentially replacing plot_gaia_cmd.py and potentially 
displacing plot_data_checks.py
"""


from __future__ import print_function
import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, QTable
import matplotlib.pyplot as plt
import scipy.stats as scistats

import plotting_dicts as pod


plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
#plt.rc('font', size = 11)
plt.rc('lines', markersize = 5)


bcolor='g'
ncolor= 'cyan'
list_color = '#1ca1f2'
annotate= True #controls whether or not object names appear beside points in the scatter plots. Should be turned off for >~20 targets appearing close together
target_label=''
mc_number = 10000
percent_off = 34 #1-sigma equivalent


zeropoint_dict={"g": [25.6883657251, 0.0017850023],
                "bp": [ 25.3513881707 , 0.0013918258],
                "rp": [24.7619199882, 0.0019145719]} #from Evans et al 2018, the DR2 values [ZP, sigma]

#############

#target_input='20190829_alkaliWD_targeted_gaia_scbd.csv'
target_input='20190829_DZNas.csv'

#####

other_target_input='20190516B_retargeted_purple_search_gaia_scbd.csv'


###########################

target_table = Table.read(target_input)
other_target_table=Table.read(other_target_input)



def get_mag(flux, filter_string):
    mag0 = zeropoint_dict[filter_string][0]
    return -2.5*np.log10(flux) +mag0

def get_mc_distribution(value, error):
    error_distribution = np.random.normal(loc= value, scale = error, size = mc_number)
    return error_distribution


def remove_negative(array, verbose= True):
    output_array = array[np.where(array>0)]
    if (verbose and array.shape[0]-output_array.shape[0] >0):
        print('Removed ' +str(array.shape[0]-output_array.shape[0]) + ' negatives')
    return output_array

def match_sizes(change_array, match_array):
    """
    Intended to keep compatibility with an array that has had negatives removed
    """
    try:
        min_inds = np.nanmin([change_array.shape[0], match_array.shape[0]])
        
        return change_array[:min_inds], match_array[:min_inds]
    except AttributeError:
        #the inputs aren't actually arrays
        return change_array, match_array
    
    
def get_errors(distribution, percent_off = percent_off):
    """
    values for the error bars on the plot
    
    Returns
    [- value, + value], so to get the points on the plot where they'd be located you do
    value - [- value] , value + [- value]
    Basically these are the width of the uncertainty range on either side.
    """
    low_bar = np.nanpercentile(distribution, 50-percent_off)
    median = np.nanmedian(distribution)
    high_bar = np.nanpercentile(distribution, 50+percent_off)
    try:
        return np.array([[median-low_bar],[high_bar-median]])
    #except astropy.units.core.UnitsError as error:
    except u.core.UnitsError as error:
        return  np.array([[median.value-low_bar],[high_bar-median.value]])
    

def get_bp_rp_excess(row, color= list_color):
    bp_string= 'phot_bp_mean_flux'
    rp_string= 'phot_rp_mean_flux'
    g_string= 'phot_g_mean_flux'
    calc_excess= (row[bp_string]+row[rp_string])/row[g_string]
    rp_dist= get_mc_distribution(row[rp_string], row[rp_string+'_error'])
    bp_dist= get_mc_distribution(row[bp_string], row[bp_string+'_error'])
    g_dist= get_mc_distribution(row[g_string], row[g_string+'_error'])
    calc_excess_dist= (bp_dist+rp_dist)/g_dist
    excess_error= get_errors(calc_excess_dist)
    #plt.hist(calc_excess_dist, bins =50)
    #plt.title(row['name'])
    #plt.xlabel('bp rp excess')
    #plt.show()
    rp_mag= get_mag(row[rp_string], 'rp')
    rp_mag_dist= get_mag(rp_dist, 'rp')
    bp_mag= get_mag(row[bp_string], 'bp')
    bp_mag_dist= get_mag(bp_dist, 'bp')
    bp_rp= bp_mag-rp_mag
    bp_rp_dist= bp_mag_dist-rp_mag_dist
    bp_rp_error= get_errors(bp_rp_dist)
    
    plt.errorbar(np.copy(bp_rp), np.copy(calc_excess), yerr = np.copy(excess_error),  xerr = np.copy(bp_rp_error), marker = 'o', markersize = 6, color = color, capsize = 4, label = target_label, linestyle ='none')
    print('excess_error.shape', excess_error.shape)
    print('bp_rp_error.shape', bp_rp_error.shape)
    if annotate:
            #plt.annotate(str(row['name']),xy=(bp_rp+bp_rp_error[1,0], calc_excess+excess_error[1,0]), xycoords='data', xytext=(bp_rp+bp_rp_error[1,0], calc_excess+excess_error[1,0]), textcoords= 'data' , fontsize=8, color =color)
            plt.annotate(str(row['name']),xy=(bp_rp ,calc_excess), xycoords='data', xytext=(bp_rp, calc_excess), textcoords= 'data' , fontsize=8, color =color)
    #rp_error= get_errors(rp_dist)
    #bp_error=get_errors(bp_dist)
    #g_error= get_errors(g_dist)
    
    return


def plot_target_table(input_table, x_variable, y_variable, color=list_color):
    for row in input_table:
        get_bp_rp_excess(row, color=color)
    plt.xlabel('BP-RP')
    plt.ylabel('phot_bp_rp_excess_factor')
    plt.legend()
    plt.title(target_input)
    #plt.show()
    
    
    return


plot_target_table(other_target_table, 'a', 'b', color='r')
plot_target_table(target_table, 'a', 'b')
xvals= np.linspace(-2,6.,1000)
yvals= 1.65-0.03*(xvals-2.2)**2. +0.1*xvals
y2vals= 1.7+0.06*xvals**2.
yvalsL= 1.3+0.06*(xvals)**2.
yvalsL2=1.0+0.015*xvals**2.
plt.plot(xvals,yvalsL, linestyle='--', color='magenta', label = "Lindegren's cut")
plt.plot(xvals,yvalsL2,linestyle='--',color='magenta')
#plt.plot(xvals,yvals, linestyle='--', color='g', label= "Ben's cut")
plt.plot(xvals, y2vals, linestyle = '--', color='cyan', label = "Nicola's cut")
plt.axhline(y=1.8,linestyle='--', color='blue', label='1.8 flat cut')
plt.legend(loc='best')
plt.show()



