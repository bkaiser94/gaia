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

#cmap='viridis'




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
grid_num = 225

#basis vector components for orthnormal bp vs. rp plot
v1bp=-0.78013779
v1rp=-0.625607727

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
    #generic_input= 'Lindegren_appC_selA_obsnum.csv'
    #generic_input = 'Lindegren_washy_cloud_bigger.csv'
    #generic_input='Lindegren_appC_selB_washy_cloud_bigger.csv'
    #generic_input= 'Lindegren_appC_selB_antiC_nobulgedisk.csv'
    #generic_input= 'Lindegren_appC_selB_antiC_cut2.csv'
    #generic_input='Lindegren_appC_selB_antiC_cut2_gaia_sc.csv'
    #generic_input='20190121_excess_interesting_gaia.csv'
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
#########################################
col_pairs=[
    ['bp_rp','mg'],
    ['astrometric_pseudo_colour','mg'],
    ['ra','dec'],
    ['l','b'],
    ['parallax','phot_g_mean_mag'],
    ['phot_bp_mean_mag','phot_rp_mean_mag'],
    ['bp_rp','g_rp'],
    ['bp_rp', 'bp_g'],
    ['bp_rp','astrometric_pseudo_colour'],
    ['astrometric_excess_noise','phot_g_mean_mag'],
    ['astrometric_excess_noise','phot_bp_mean_mag'],
    ['astrometric_excess_noise','phot_rp_mean_mag'],
    ['astrometric_excess_noise','phot_bp_rp_excess_factor'],
    ['bp_rp', 'phot_bp_rp_excess_factor'],
    ['phot_g_mean_mag', 'phot_bp_rp_excess_factor'],
    ['phot_bp_mean_mag', 'phot_bp_rp_excess_factor'],
    ['phot_rp_mean_mag', 'phot_bp_rp_excess_factor'],
    ['ra','pmra'],
    ['dec','pmra'],
    ['dec','pmdec'],
    ['ra','pmdec'],
    ['pmra','pmdec']]


############################################
try:
    cprime_y= v1bp * generic_table['phot_bp_mean_mag'] + v1rp*generic_table['phot_rp_mean_mag']
    cprime_x= v1rp * generic_table['phot_bp_mean_mag'] - v1bp*generic_table['phot_rp_mean_mag']
except KeyError as error:
    print("Don't have", error, '\nSo no alternative basis colors')
    
    
def scatter_plot(string_pair):
    if string_pair[0]=='parallax':
        x_array=1000./generic_table['parallax']
        string_pair[0]= '1000/parallax'
        y_array = generic_table[string_pair[1]]
    elif string_pair[1]=='parallax':
        y_array=1000./generic_table['parallax']
        string_pair[1]= '1000/parallax'
        x_array = generic_table[string_pair[0]]
    if string_pair[0]=='astrometric_pseudo_colour':
        x_array=1./generic_table['astrometric_pseudo_colour']
        string_pair[0]= '1/astrometric_pseudo_colour'
        y_array = generic_table[string_pair[1]]
    else:
        x_array = generic_table[string_pair[0]]
        y_array = generic_table[string_pair[1]]
    try:
        #plt.scatter(x_array, y_array, c= generic_table['phot_proc_mode'], alpha=0.5)
        sorted_order= np.argsort(generic_table['phot_bp_rp_excess_factor'])
        sorted_x_array=x_array[sorted_order]
        sorted_y_array= y_array[sorted_order]
        #plt.scatter(x_array, y_array, c= generic_table['phot_bp_rp_excess_factor'], alpha=0.5, markersize=4)
        plt.scatter(sorted_x_array, sorted_y_array, c= generic_table['phot_bp_rp_excess_factor'][sorted_order], alpha=1, s=8, edgecolor='none')
    except KeyError as error:
        print('No ', error,'\nTherefore using uniform color for scatter plot.')
        plt.scatter(x_array, y_array, alpha=0.5)
    #plt.plot([5.33,22.3342],[5.49,19.126])
    plt.xlabel(string_pair[0])
    plt.ylabel(string_pair[1])
    plt.title(string_pair[1] + ' vs. ' + string_pair[0])
    plt.show()
    return

def hexbin_plot(string_pair):
    if string_pair[0]=='parallax':
        x_array=1000./generic_table['parallax']
        string_pair[0]= '1000/parallax'
        y_array = generic_table[string_pair[1]]
    elif string_pair[1]=='parallax':
        y_array=1000./generic_table['parallax']
        string_pair[1]= '1000/parallax'
        x_array = generic_table[string_pair[0]]
    else:
        x_array = generic_table[string_pair[0]]
        y_array = generic_table[string_pair[1]]
    polything = plt.hexbin(x_array,y_array, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1)
    counts = polything.get_array()
    #counts= np.sqrt(counts)
    counts=np.log(counts)
    polything.set_array(counts)
    polything.autoscale()
    plt.xlabel(string_pair[0])
    plt.ylabel(string_pair[1])
    plt.title(string_pair[1] + ' vs. ' + string_pair[0])
    #plt.plot([5.33,22.3342],[5.49,19.126])
    plt.show()
    return

def bp_rp_cut_line(string_pair):
    if ((string_pair[0]== 'bp_rp') and (string_pair[1]=='phot_bp_rp_excess_factor')):
        xvals= np.linspace(-2,6.,1000)
        yvals= 1.3+0.06*(xvals)**2.
        plt.plot(xvals,yvals, linestyle='--', color='magenta')
    else:
        pass
    return


plt.scatter(cprime_x, cprime_y, alpha=0.5)
#polything = plt.hexbin(cprime_x,cprime_y, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1)
#counts = polything.get_array()
#counts= np.sqrt(counts)
##counts=np.log(counts)
#polything.set_array(counts)
#polything.autoscale()
plt.xlabel('cprime_x')
plt.ylabel('cprime_y')
plt.show()

plt.scatter(cprime_x, generic_table['phot_bp_rp_excess_factor'], alpha=0.5)
#polything = plt.hexbin(cprime_x,cprime_y, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1)
#counts = polything.get_array()
#counts= np.sqrt(counts)
##counts=np.log(counts)
#polything.set_array(counts)
#polything.autoscale()
plt.xlabel('cprime_x')
plt.ylabel('phot_bp_rp_excess_factor')
plt.show()

#scatter_plot(['phot_bp_mean_mag','phot_rp_mean_mag'])
#hexbin_plot(['phot_bp_mean_mag','phot_rp_mean_mag'])

    
for string_pair in col_pairs:
    if string_pair[1]=='mg':
        plt.gca().invert_yaxis()
    else:
        pass
    try:
        bp_rp_cut_line(string_pair)
        scatter_plot(string_pair)
        #hexbin_plot(string_pair)
    except KeyError as error:
        print('No column named', error, '\nSkipping', string_pair[1] + ' vs. ' + string_pair[0])
       




plt.hist(generic_table['phot_proc_mode'])
plt.title('phot_proc_mode')
plt.show()

plt.hist(generic_table['astrometric_excess_noise'], bins= 100)
plt.title('astrometric_excess_noise')
plt.show()

