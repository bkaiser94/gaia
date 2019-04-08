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
#########3
ncolor= 'cyan' #Nicola's line color
bcolor= 'g' #Bens' line color
lcolor= 'magenta' #Lindegren's line color


#######3error distribution variables
mc_number = 10000
percent_off = 34 #1-sigma equivalent
#############


#num_targs = 'all'
#num_targs = '47Tuc'
num_targs= 'Lindegren'
selection_letter= 'B'
distance = 200
grid_num = 100

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
    #generic_input   ='usdMs_gaia.csv'
    #generic_input= 'Lindegren_odd_survivors_gaia_sc.csv'
    #generic_input= 'massive_zzceti_gaia.csv'

    #generic_input='NaD_objects.csv'
    generic_input='WD_cooling_tip_gaia.csv'

    #generic_input='exc1_8_2_2_purple_search_gmaglimit_gaia_sc.csv'
    #generic_input='Lindegren_appB_bulge_only.csv'
    #generic_input='20190107_chris_merge_gaia.csv'
    #generic_input= '20190128_wdMS_gaia.csv'
    #generic_input= 'RNe_gaia.csv'
    #generic_input= 'coolDZ_Na_gaia.csv'
    #generic_input= 'Lindegren_appC_selA_obsnum.csv'
    #generic_input = 'Lindegren_washy_cloud_bigger.csv'
    #generic_input='Lindegren_appC_selB_washy_cloud_bigger.csv'
    #generic_input= 'Lindegren_appC_selB_antiC_nobulgedisk.csv'
    #generic_input= 'Lindegren_appC_selB_antiC_cut2.csv'
    #generic_input='Lindegren_appC_selB_antiC_cut2_gaia_sc.csv'
    #generic_input='20190121_excess_interesting_gaia.csv'
    #generic_input= 'Lindegren_appC_selA_hv_region.csv'
    #generic_input='Lindegren_appC_altC_noBDLMC.csv'
    #generic_input= '20190107_chris_merge_gaia.csv'
    #generic_input= 'ar_sco_gaia.csv'
    #generic_input= '20190123_new_red_things_gaia_sc.csv'
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
#########################################
col_pairs=[
    ['bp_rp','mg'],
    ['g_rp','mg'],
    ['astrometric_pseudo_colour','mg'],
    ['bp_rp', 'phot_bp_rp_excess_factor'],
    ['mg','mean_varpi_factor_al'],
    ['phot_bp_rp_excess_factor','mean_varpi_factor_al'],
    ['astrometric_excess_noise','mean_varpi_factor_al'],
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
    ['phot_g_mean_mag', 'phot_bp_rp_excess_factor'],
    ['phot_bp_mean_mag', 'phot_bp_rp_excess_factor'],
    ['phot_rp_mean_mag', 'phot_bp_rp_excess_factor'],
    ['ra','pmra'],
    ['dec','pmra'],
    ['dec','pmdec'],
    ['ra','pmdec'],
    ['pmra','pmdec']]

col_singles=[
    'ra',
    'dec',
    'phot_bp_rp_excess_factor',
    'phot_proc_mode',
    'astrometric_excess_noise',
    'parallax']
    

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
    elif string_pair[0]=='astrometric_pseudo_colour':
        x_array=1./generic_table['astrometric_pseudo_colour'] * 1e4
        string_pair[0]= '1/astrometric_pseudo_colour (angstroms)'
        y_array = generic_table[string_pair[1]]
    elif string_pair[1]=='astrometric_pseudo_colour':
        y_array=1./generic_table['astrometric_pseudo_colour'] * 1e4
        string_pair[1]= '1/astrometric_pseudo_colour (angstroms)'
        x_array = generic_table[string_pair[0]]
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
    #plt.title(string_pair[1] + ' vs. ' + string_pair[0])
    plt.title(generic_input)
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
    elif string_pair[0]=='astrometric_pseudo_colour':
        x_array=1./generic_table['astrometric_pseudo_colour'] * 1e4
        mask= np.where(np.abs(x_array)>20000)
        x_array[mask]=0
        string_pair[0]= '1/astrometric_pseudo_colour (angstroms)'
        y_array = generic_table[string_pair[1]]
    elif string_pair[1]=='astrometric_pseudo_colour':
        y_array=1./generic_table['astrometric_pseudo_colour'] * 1e4
        mask= np.where(np.abs(y_array)>20000)
        y_array[mask]=0
        string_pair[1]= '1/astrometric_pseudo_colour (angstroms)'
        x_array = generic_table[string_pair[0]]
    else:
        x_array = generic_table[string_pair[0]]
        y_array = generic_table[string_pair[1]]
        #y_array =np.log10(y_array)
    polything = plt.hexbin(x_array,y_array, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1)
    counts = polything.get_array()
    #counts= np.sqrt(counts)
    counts=np.log(counts)
    polything.set_array(counts)
    polything.autoscale()
    plt.xlabel(string_pair[0])
    plt.ylabel(string_pair[1])
    #plt.title(string_pair[1] + ' vs. ' + string_pair[0])
    plt.title(generic_input)
    #plt.plot([5.33,22.3342],[5.49,19.126])
    plt.show()
    return

def hist_plot(string_single):
    plt.hist(generic_table[string_single])
    plt.title(string_single)
    plt.show()
#def bp_rp_cut_line(string_pair):
    #if ((string_pair[0]== 'bp_rp') and (string_pair[1]=='phot_bp_rp_excess_factor')):
        #xvals= np.linspace(-2,6.,1000)
        #yvals= 1.3+0.06*(xvals)**2.
        #yvals2=1.0+0.015*xvals**2.
        #plt.plot(xvals,yvals, linestyle='--', color='magenta', label = "Lindegren's cut")
        #plt.plot(xvals,yvals2,linestyle='--',color='magenta')
    #else:
        #pass
    #return


def plot_cut_lines(string_pair):
    #colour excess cuts
    if ((string_pair[0]== 'bp_rp') and (string_pair[1]=='phot_bp_rp_excess_factor')):
        xvals= np.linspace(-2,6.,1000)
        yvals= 1.65-0.03*(xvals-2.2)**2. +0.1*xvals
        y2vals= 1.7+0.06*xvals**2.
        yvalsL= 1.3+0.06*(xvals)**2.
        yvalsL2=1.0+0.015*xvals**2.
        plt.plot(xvals,yvalsL, linestyle='--', color='magenta', label = "Lindegren's cut")
        plt.plot(xvals,yvalsL2,linestyle='--',color='magenta')
        plt.plot(xvals,yvals, linestyle='--', color='g', label= "Ben's cut")
        plt.plot(xvals, y2vals, linestyle = '--', color='cyan', label = "Nicola's cut")
        plt.legend()
    if ((string_pair[0]== 'bp_rp') and (string_pair[1]=='mg')):
        x1vals= np.linspace(-1,-0.184268, 300)
        y1vals= np.ones(x1vals.shape)*5
        x2vals= np.linspace(np.max(x1vals),0.297505,300)
        y2vals= 5.93+5.047*x2vals
        x3vals= np.linspace(np.max(x2vals), 1.7, 300)
        y3vals = 6*x3vals**3.-21.77*x3vals**2.+27.91*x3vals+0.897
        y4vals= np.linspace(14.9067,16,300)
        x4vals= np.ones(y4vals.shape)*1.7
        plt.plot(x1vals,y1vals, color= ncolor, label="Nicola's cut")
        plt.plot(x2vals, y2vals, color= ncolor)
        plt.plot(x3vals, y3vals, color= ncolor)
        plt.plot(x4vals, y4vals, color= ncolor)
        plt.legend()
    if  (string_pair[1]=='mean_varpi_factor_al'):
        plt.axhline(y=-0.23, color='magenta', label='acceptable range', linestyle='--')
        plt.axhline(y=0.32,color='magenta', linestyle='--')
        plt.axhline(y=-0.73, color='k', label='"constrained range"')
        plt.axhline(y=0.73, color='k')
        plt.legend()
    else:
        pass
    return

def pseudo_colour_bp_rp_line(string_pair):
    if ((string_pair[0] == 'bp_rp') and (string_pair[1] == 'astrometric_pseudo_colour')):
        xvals= np.linspace(-2,6.,1000)
        yvals = 2.0-1.8/np.pi * np.arctan(0.331+0.572*xvals-0.014*xvals**2+0.045*xvals**3)
        yvals= 1/yvals * 1e4
        plt.plot(xvals,yvals,color='magenta')
    else:
        pass
    return



#plt.scatter(cprime_x, cprime_y, alpha=0.5)
polything = plt.hexbin(cprime_x,cprime_y, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1)
counts = polything.get_array()
#counts= np.sqrt(counts)
counts=np.log(counts)
polything.set_array(counts)
polything.autoscale()
plt.xlabel('cprime_x')
plt.ylabel('cprime_y')
plt.show()

#plt.scatter(cprime_x, generic_table['phot_bp_rp_excess_factor'], alpha=0.5)
polything = plt.hexbin(cprime_x,generic_table['phot_bp_rp_excess_factor'], gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1)
counts = polything.get_array()
#counts= np.sqrt(counts)
counts=np.log(counts)
polything.set_array(counts)
polything.autoscale()
plt.xlabel('cprime_x')
plt.ylabel('phot_bp_rp_excess_factor')
plt.show()


try:
    calc_pseudo = 2.0-1.8/np.pi * np.arctan(0.331+0.572*generic_table['bp_rp']-0.014*generic_table['bp_rp']**2+0.045*generic_table['bp_rp']**3)
    calc_pseudo= 1/calc_pseudo * 1e4
    #polything = plt.hexbin(calc_pseudo-1/generic_table['astrometric_pseudo_colour']*1e4,generic_table['mg'], gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1)
    #counts = polything.get_array()
    polything = plt.scatter(calc_pseudo-1/generic_table['astrometric_pseudo_colour']*1e4,generic_table['mg'], s=4, alpha=0.1, edgecolor='none')
    #counts= np.sqrt(counts)
    #counts=np.log(counts)
    #polything.set_array(counts)
    #polything.autoscale()
    plt.xlabel(r'$\lambda_{eff, bp-rp} - \lambda_{eff, pseudo}$')
    plt.gca().invert_yaxis()
    plt.ylabel('mg')
    plt.show()
except KeyError:
    pass

#scatter_plot(['phot_bp_mean_mag','phot_rp_mean_mag'])
#hexbin_plot(['phot_bp_mean_mag','phot_rp_mean_mag'])

    
for string_pair in col_pairs:
    if string_pair[1]=='mg':
        plt.gca().invert_yaxis()
    else:
        pass
    try:
        #bp_rp_cut_line(string_pair)
        plot_cut_lines(string_pair)
        pseudo_colour_bp_rp_line(string_pair)
        scatter_plot(string_pair)
        #hexbin_plot(string_pair)
    except KeyError as error:
        plt.clf()
        print('No column named', error, '\nSkipping', string_pair[1] + ' vs. ' + string_pair[0])
       

for single_string in col_singles:
    try:
        #bp_rp_cut_line(string_pair)
        hist_plot(single_string)
    except KeyError as error:
        plt.clf()
        print('No column named', error, '\nSkipping', col_singles, ' histogram')


plt.hist(generic_table['phot_proc_mode'])
plt.title('phot_proc_mode')
plt.show()

plt.hist(generic_table['astrometric_excess_noise'], bins= 100)
plt.title('astrometric_excess_noise')
plt.show()

try:
    plt.hist(generic_table['mean_varpi_factor_al'], bins= 100)
    plt.title('mean_varpi_factor_al')
    plt.axvline(x=-0.23, color='magenta', label='acceptable range', linestyle='--')
    plt.axvline(x=0.32,color='magenta', linestyle='--')
    plt.axvline(x=-0.73, color='k', label='"constrained range"')
    plt.axvline(x=0.73, color='k')
    mean_val= np.nanmean(generic_table['mean_varpi_factor_al'])
    med_val=np.nanmedian(generic_table['mean_varpi_factor_al'])
    per1= np.percentile(generic_table['mean_varpi_factor_al'], 1)
    print('1st percentile: ', per1)
    per99= np.percentile(generic_table['mean_varpi_factor_al'], 99)
    print('99th percentile: ', per99)
    plt.axvline(x=mean_val, label='mean: ' + str(mean_val)[:5], linestyle='--', color='cyan')
    plt.axvline(x=med_val, label='med: '+ str(med_val)[:5], linestyle='--', color='red')
    plt.legend()
    plt.show()
except KeyError:
    pass

