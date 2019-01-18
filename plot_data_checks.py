"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-01-17

Plot various different gaia data columns against each other or as histograms to check for correlations or patterns
in a given sample of gaia data. Also allows for checking of the phot_proc_mode, if it's included in the table that is
used for the plots...

"""



from __future__ import print_function
import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, QTable
import matplotlib.pyplot as plt
import scipy.stats as scistats
#import seaborn as sns
import astropy


#import passband_model_convolution as pmc
import gaia_extinction
#import wdatmos
import plotting_dicts as pod


plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
#plt.rc('font', size = 11)
plt.rc('lines', markersize = 5)

absmag_band= 'g'
colours= ['g','rp']





#axes_x= [-1.5, 6]
#axes_y = [-4,18]


axes_x= [-1, 6]
axes_y = [-5,25]

y_fill= axes_y[1]
x_fill= axes_x[0]

list_color = '#1ca1f2'
single_list=True#turns off the original for-loop method for plotting from a single list
error_bar=True #turns off error bars on the multiple list plot, meaning it has no effect on anything if single_list==True
annotate= False #controls whether or not object names appear beside points in the scatter plots. Should be turned off for >~20 targets appearing close together
parallax_correction = 0.029 #from Lindgren et al 2018


#######3error distribution variables
mc_number = 10000
percent_off = 34 #1-sigma equivalent
#############

target_input='20190107_chris_merge_gaia.csv'
target_input='l-0.3bp_g_gaia_corr_full.csv'
target_input='elm_survey_gaia.csv'
target_input='pre_elms_gaia.csv'
target_input = '20190109_blue_gaia.csv'
target_input='hot_wind_wds_gaia.csv'
target_input='Eriks_disk_candidates_gaia.csv'
target_input='dC_sample_roulston2018_gaia.csv'
target_input='20190111_red_things_gaia.csv'
target_input='alt_red_things_g_rp_greater_17_gaia.csv'
target_input='WD_cooling_tip_gaia.csv'
target_input='weird_CPM_binary_gaia.csv'
target_label= ''

num_targs = 'all'
#num_targs = '47Tuc'
num_targs= 'Lindegren'
selection_letter= 'C'
distance = 200
grid_num = 500


#################################################################3
########## End of things that should be edited for a given run######################
#######################################################################
#Selecting the csv file to use to generate the background CMD
if num_targs == 'all':
    print('Distance-limited Sample like Figure 6 from DR2HRD')
    generic_input = 'all_'+str(distance)+'pc_gaia_corr.csv'
    title_suffix = ' in the ' + str(distance)+ 'pc Gaia DR2 CMD'
elif num_targs== '47Tuc':
    generic_input= "47Tuc_10arcmin.csv"
elif num_targs== 'Lindegren':
    generic_input = 'Lindegren_appC_sel'+selection_letter + '.csv'
    generic_input= 'Lindegren_appC_selA_obsnum.csv'
    generic_input = 'Lindegren_washy_cloud_bigger.csv'
    generic_input='Lindegren_appC_selB_washy_cloud_bigger.csv'
    generic_input= 'Lindegren_appC_selB_antiC_nobulgedisk.csv'
    #generic_input= 'Lindegren_appC_selA_hv_region.csv'
else:
    num_targs = int(num_targs)
    generic_input = 'top'+str(num_targs) + '_nearby_gaia.csv'
    title_suffix = 'in the ' +str(num_targs)+ ' star sample following DR2HRD Figure 1'
    
############################

zeropoint_dict={"g": [25.6883657251, 0.0017850023],
                "bp": [ 25.3513881707 , 0.0013918258],
                "rp": [24.7619199882, 0.0019145719]} #from Evans et al 2018, the DR2 values [ZP, sigma]


################################
#Reading in the tables for the background and target files
generic_table = Table.read(generic_input)
target_table = Table.read(target_input)

col_pairs=[
    ['parallax','phot_g_mean_mag'],
    ['phot_bp_mean_mag','phot_rp_mean_mag'],
    ['astrometric_excess_noise','phot_g_mean_mag'],
    ['bp_rp','g_rp']]

plt.scatter(1000./generic_table['parallax'], generic_table['phot_g_mean_mag'])
plt.show()
plt.scatter(generic_table['phot_bp_mean_mag'], generic_table['phot_rp_mean_mag'])
plt.show()

#plt.scatter(generic_table['phot_bp_mean_mag'], generic_table['phot_rp_mean_mag'])
#plt.show()


plt.scatter(generic_table['astrometric_excess_noise'], generic_table['phot_g_mean_mag'])
plt.title(' phot_g_mean_mag vs. astrometric_excess_noise')
plt.show()

plt.scatter(generic_table['bp_rp'], generic_table['g_rp'])
plt.title(' bp_rp vs. g_rp')
plt.show()

plt.scatter(generic_table['astrometric_excess_noise'], generic_table['phot_bp_mean_mag'])
plt.title(' phot_bp_mean_mag vs. astrometric_excess_noise')
plt.show()

plt.scatter(generic_table['astrometric_excess_noise'], generic_table['phot_rp_mean_mag'])
plt.title(' phot_rp_mean_mag vs. astrometric_excess_noise')
plt.show()


plt.hist(generic_table['phot_proc_mode'])
plt.title('phot_proc_mode')
plt.show()

plt.hist(generic_table['astrometric_excess_noise'], bins= 100)
plt.title('astrometric_excess_noise')
plt.show()

plt.scatter(generic_table['ra'],generic_table['dec'])
plt.title('ra vs. dec')
plt.show()


plt.scatter(generic_table['l'],generic_table['b'])
plt.title('l vs. b')
plt.show()
