"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-01-07 (spun-off from plot_alt_cmd.py on 2022-11-16)

This script should plot lists of targets over the CMD for a given distance, and allow for the creation of different
CMDs than M_G vs. BP-RP, i.e. M_G vs. G-RP

It's going to borrow heavily from plot_gaia_cmd.py, but I'm going to try to keep it clean of all of the extra stuff
I have in there to do calculations of a single target.

Ok. This new (2022-11-16) script should be able to handle the el-Badry files with 2 sets of targets to plot simultaneously, so I should probably just add a kwarg to provide the number as string to append to whatever col names in the table for the existing functions.

I also need to create a background CMD set to pull from for the newer data releases because 
the bandpasses all changed (allegedly), so presumably the CMD should have different shapes 
and values.

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
#import astropy


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





axes_x= [-1.5, 6]
#axes_x= [-0.5, 2]
axes_y = [-4,18]


#axes_x= [-3, 6]
#axes_y = [-5,25]

y_fill= axes_y[1]
x_fill= axes_x[0]

list_color = '#1ca1f2'
single_list=False#turns off the original for-loop method for plotting from a single list
error_bar=True #turns off error bars on the multiple list plot, meaning it has no effect on anything if single_list==True
annotate= True #controls whether or not object names appear beside points in the scatter plots. Should be turned off for >~20 targets appearing close together
#parallax_correction = 0.029 #from Lindgren et al 2018
parallax_correction=0.

#bcolor='g'
#ncolor= 'cyan'

markersize=4

bcolor='purple'
ncolor= 'red'

#######3error distribution variables
mc_number = 10000
percent_off = 34. #1-sigma equivalent
#percent_off = 99.7/2. #3-sigma equivalent
#############

#target_input='20190107_chris_merge_gaia.csv'
#target_input= '20190516_targeted_purple_search_gaia_scnb.csv'
#target_input= 'NLTT5306_comps_gaia.csv'
#target_input='20190511_observed_objects_gaia.csv'
#target_input= '20190405_purple_search_gaia_sc.csv'
#target_input= 'mdwarf_spTs_gaia_sc.csv'
#target_input= 'Lindegren_appC_selC_noBLMC_wd1401_bincomp_gaia_sc.csv'
#target_input= 'Lindegren_appC_selC_noBLMCbperr_wd1401_bincomp_gaia_sc.csv'
#target_input='SDSSJ0744p4649_gaia.csv'
#target_input= 'WDJ0205m053ultracool_gaia.csv'
#target_input='top100_blue_cmd.csv'
#target_input='20190109_star1_ultracool.csv'
#target_input='20210305B_ultracool_switchback_gaia_gaia_scbd.csv'
#target_input= 'observed_purple_400m2.csv'
#target_input= 'exc1_8_2_2_purple_search_gmaglimit_gaia_sc.csv'
#target_input= 'expanded_purple_search_gmaglimit_gaia_sc.csv'
#target_input='20190422_obs_objects_gaia.csv'
#target_input= 'sdssj1330p6435_gaia.csv'
#target_input= 'sdssj1330_similar_gaia_sc.csv'
#target_input= 'sdssj1330_similar_subset_gaia_sc.csv'
#target_input= 'sdssj1330_similar_subset_observed.csv'
#target_input ='20190516B_retargeted_purple_search_gaia_scbd.csv'
#target_input= '20190730_obs_objects.csv'
#target_input='GaiaJ1453m2258_andnearones.csv'
#target_input= 'gaiaj1644m0449_gaia.csv'
#target_input='20190829_alkaliWD_targeted_gaia_scbd.csv'
#target_input='20190829_DZNas.csv'
#target_input= 'BU1941m1919_gaia.csv'
#target_input='BU1941m1919_gaia_edr3.csv'
#target_input='20210305_DZs_wBe_and_J1113_gaia.csv'
#target_input='WDJ0850p1956_gaia_LiWDcand.csv'
#target_input='HD113083_omegacen_gaiadr2.csv'
#target_input='PSO_ucool_gaia.csv'
#target_input='WD_DZLi.csv'
#target_input='20220929_reallydim_WD_search_gaia_scbd.csv'
#target_input= 'SDSSJ1029p1729_extlowmetal_gaia.csv'
#target_input='20190616_TIC294.csv'
#target_input= '20190901and02_obs_objects.csv'
#target_input='20190616_TIC294andother.csv'
#target_input= 'josh_object.csv'
#target_input= '20190601_obs_objects.csv'
#target_input= '20190528_named_objects_from_retarg.csv'
#target_input='sdssj1240p6710_DS_oxygenrich.csv'
#target_input= 'psrj1048p2339_gaia.csv'
#target_input='20190516B_retargeted_purple_search_gaia_scbd_20220804_update.csv'
#target_input='BPSCS29528m0028_CEMPs_gaiaDR2.csv'
#target_input='CEMPs_Narich_gaiaDR2.csv'
#target_input='dimWDMS_elBadry_eDR3.fits'
#target_input='dimWDMS_F0toK7_eDR3.fits'
#target_input='dimWDMS_F0toK7_eDR3_highconf.fits'
#target_input='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf.csv'
#target_input='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf_Vincent2023bClass.csv'
#target_input='2MASSJ0916m4215_DZLi_Vennes_gaiadr2.csv'
#target_input='2MASSJ0916m4215_DZLi_Vennes_gaiaedr3.csv'
#target_input='WDJ0212m5522_coolDZ_gaiaDR3.csv'
#target_input='sdss_kstars_fehlm05_xmatch_gaiadr3.csv'
#target_input='WDs_mistaken_for_Kstars_by_SDSS_gaiaDR3names.csv'
#target_input='WDs_mistaken_for_Kstars_nofehlimit_by_SDSS_gaiaDR3.csv'
#target_input='sdss_kstars_nofehlimit_xmatch_gaiadr3.csv'
target_input='SDSSJ0804p5130_gaiadr3.csv'
#target_input='20240116_1105_K5FeHLm05_15arcsecradgaiaDR3.csv'
#target_input='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf_Vincent2023bClass_20240125update.csv'
#target_input='WDJ1948m1011_gaiaDR3.csv'
#target_input='WDJ1515p1911_gaiaDR3.csv'


#other_target_input='20220929_reallydim_WD_search_gaia_scbd.csv'

other_target_input= 'mdwarf_spTs_gaia_sc.csv'
#other_target_input='WDJ1203m0012_gaiaDR3.csv'
#other_target_input='20210201_DZs_for_J1636paper_gaia_gaiaDR3.csv'
#other_target_input='gaiaeDR3_spectral_type_objects.csv'
#other_target_input='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf.csv'
#other_target_input='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf_Vincent2023bClass.csv'

#other_target_input='DESJ2147m4035_Appsweirdcool28pc_dr2.csv'
#other_target_input='20190516B_retargeted_purple_search_gaia_scbd_20220804_update.csv'
#other_target_input= 'SDSSJ1150p2403_gaia.csv'
#other_target_input= 'HE1327m2326_extlowmetal_gaia.csv'
#other_target_input= '20190730_obs_objects.csv'
#other_target_input='gaiaj1644m0449other_gaia.csv'
#other_target_input='GaiaJ1453m2258_andnearones.csv'
#other_target_input= '20190703_obs_objects.csv'
#other_target_input='sdssj1330p6435_gaia.csv'
#other_target_input='20190516B_retargeted_purple_search_gaia_scbd.csv'
#other_target_input='wd2251m070_gaia.csv'
#other_target_input= 'sdssj1330_similar_subset_observed.csv'
#other_target_input='20190616_obs.csv'

target_label= ''

num_targs = 'all'
#num_targs = '47Tuc'
#num_targs= 'Lindegren'
selection_letter= 'C'
#distance = 100
distance=200
grid_num = 225


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
    #generic_input='rand_astrometric_primaries1.csv'
    #generic_input = 'Lindegren_appC_sel'+selection_letter + '_200pc.csv'
    #generic_input='Lindegren_appC_selB_antiC_cut2_gaia_sc.csv'
    #generic_input='Lindegren_appC_altC_noBDLMC.csv'
    #generic_input= 'Lindegren_appC_selB_antiC_cut2.csv'
    #generic_input= 'Lindegren_appC_selB_antiC_cut2_gaia_sc.csv'
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
other_target_table=Table.read(other_target_input)
#target_table=target_table[np.where(target_table['repeat'] == 'False')]
#target_table=target_table[np.where(target_table['priority'] <= 4.5)]
#target_table=target_table[np.where(target_table['priority']> 98)]
#other_target_table=other_target_table[np.where(other_target_table['priority'] > 98)]
#blank_val=np.copy(target_table[2]['name'])
#target_table=target_table[np.where(target_table['name']!='')]
#target_table=target_table[np.where(target_table['wd_dm']> 0)]

#print('name',target_table[2]['name'])
####################################

#Added this part to be able to make plots to indicate which objects have continuous spectra and which don't.
#xp_cont_inds=np.where(target_table['has_xp_continuous']=="True")
#xp_nocont_inds=np.where(other_target_table['has_xp_continuous']=="False")

#print('\n\n==== xp_cont_inds')
#print(xp_cont_inds)
#print(xp_nocont_inds)
#print('\n\n=====+++++====')

#target_table=target_table[xp_cont_inds]
#other_target_table=other_target_table[xp_nocont_inds]






##########################

def plot_ben_cuts():
    x1vals= np.linspace(0.82, 0.97, 100)
    y1vals= 8.87*x1vals + 6.8
    x2vals= np.linspace(1.32, 1.52, 100)
    y2vals= 8.15*x2vals+3.312
    x3vals= np.linspace(0.82, 1.32,100)
    y3vals=np.ones(x3vals.shape)*14.07
    y4vals= np.linspace(15.4, 18, 100)
    x4vals= np.ones(y4vals.shape)*0.97
    y5vals= np.linspace(15.7, 18, 100)
    x5vals= np.ones(y5vals.shape)*1.52
    plt.plot(x1vals,y1vals, color= bcolor, label="Purple Object Survey")
    plt.plot(x2vals, y2vals, color= bcolor)
    plt.plot(x3vals, y3vals, color= bcolor)
    plt.plot(x4vals, y4vals, color= bcolor)
    plt.plot(x5vals, y5vals, color= bcolor)
    plt.legend()
    return

def plot_nicola_cuts():
    x1vals= np.linspace(-1,-0.184268, 300)
    y1vals= np.ones(x1vals.shape)*5
    x2vals= np.linspace(np.max(x1vals),0.297505,300)
    y2vals= 5.93+5.047*x2vals
    x3vals= np.linspace(np.max(x2vals), 1.7, 300)
    y3vals = 6*x3vals**3.-21.77*x3vals**2.+27.91*x3vals+0.897
    y4vals= np.linspace(14.9067,16,300)
    x4vals= np.ones(y4vals.shape)*1.7
    #plt.plot(x1vals,y1vals, color= ncolor, label="Nicola's cut")
    plt.plot(x1vals,y1vals, color= ncolor, label="Gentile Fusillo et al. 2019")
    plt.plot(x2vals, y2vals, color= ncolor)
    plt.plot(x3vals, y3vals, color= ncolor)
    plt.plot(x4vals, y4vals, color= ncolor)
    plt.legend()
    return



def distance_modulus(g_mag, distance):
    return g_mag - 5*np.log10(distance/10.)
    
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


def get_filter_vals(table, filter_string,num=""):
    flux_string = 'phot_'+filter_string+ '_mean_flux'
    phot_mean_flux = table[flux_string+num]
    error_string = flux_string + '_error'+num
    phot_mean_flux_error = table[error_string]
    flux_distribution = get_mc_distribution(phot_mean_flux, phot_mean_flux_error)
    return phot_mean_flux, flux_distribution



def get_colour_dif(table, plot_all = False, verbose =True, colours=['bp','rp'],num=""):
    """
    Inputs:
        table: an astropy table of Gaia values
    
        plot_all : boolean to decide if the histogram of the colours should be plotted
        
        verbose: boolean of whether or not to print things during the execution
        
        colours: two-element list that contains strings for the colour bandpasses to use.
        
        
    Outputs:
        colour_dif: the magnitude of colours[0]- magnitude of colours[1]
        
        colour_dif_error: the error bars arising from the distribution of the colour_dif values generated randomly
        
    """
    mean_flux0, dist0 = get_filter_vals(table, colours[0],num=num)
    mean_flux1  , dist1 = get_filter_vals(table, colours[1],num=num)
    mag0 = get_mag(mean_flux0, colours[0])
    mag1 = get_mag(mean_flux1, colours[1])
    if verbose:
        print(colours[0]+"_calc-bp_measured", mag0 - table['phot_' + colours[0]+ '_mean_mag'+num],num)
        print(colours[1]+ "_calc - rp_measured", mag1 - table['phot_' + colours[1]+'_mean_mag'+num],num)
    mag_dist0 = get_mag(dist0, colours[0])
    mag_dist1 = get_mag(dist1, colours[1])
    colour_dif = mag0- mag1
    colour_dif_dist= mag_dist0- mag_dist1
    colour_dif_error = get_errors(colour_dif_dist)
    if plot_all:
        plt.hist(bp_rp_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
        plt.axvline(np.nanmedian(colour_dif_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
        plt.axvline(np.nanpercentile(colour_dif_dist, 84), color = 'cyan')
        #plt.errorbar(colour_dif, 0.5, xerr = colour_dif_error, marker = '*', markersize = 8, color = 'b', label = colours[0]+"-" + colours[1], capsize = 4)
        plt.errorbar(colour_dif, 0.5, xerr = colour_dif_error, marker = '*', markersize = 8, color = 'b', label = 'changed', capsize = 4)
        #plt.xlabel(r'$G_{BP}-G_{RP}$')
        plt.xlabel(colours[0]+'-'+colours[1])
        plt.legend()
        plt.show()
    return colour_dif, colour_dif_error




def get_pass_abs_mag(table, plot_all = False, passband_string= 'g', verbose = True,num="",use_primary_parallax=False):
    mean_flux, flux_dist = get_filter_vals(table, passband_string,num)
    mag = get_mag(mean_flux, passband_string)
    flux_dist= remove_negative(flux_dist, verbose=verbose)
    mag_dist= get_mag(flux_dist, passband_string)
    if use_primary_parallax:
        parallax = table['parallax'+"1"]+parallax_correction
        parallax_error = table['parallax_error'+"1"]*1e-3
    else:
        parallax = table['parallax'+num]+parallax_correction
        parallax_error = table['parallax_error'+num]*1e-3
    #parallax = table['parallax'+num]+parallax_correction
    parallax = parallax*1e-3
    distance = 1./parallax
    #parallax_error = table['parallax_error'+num]*1e-3
    parallax_dist = get_mc_distribution(parallax, parallax_error)
    parallax_dist = remove_negative(parallax_dist, verbose= verbose)
    if verbose:
        print(passband_string+ "_calc" + "-" + passband_string+ "_measured", mag - table['phot_' +passband_string+'_mean_mag'+num],num)
    else:
        pass
    if parallax < 0:
        parallax_median = np.nanmedian(parallax_dist)
        if verbose:
            print("PARALLAX < 0!", parallax, "setting to median of positive distribution:", parallax_median)
        parallax = parallax_median
    else:
        pass
    distance = 1./parallax
    distance_dist = 1./parallax_dist
    index_length = distance_dist.shape[0]
    print("index_length",index_length)
    print("mag_dist.shape", mag_dist.shape)
    mag_dist = mag_dist[:index_length]
    mag_dist, distance_dist= match_sizes(mag_dist, distance_dist)
    abs_mag = distance_modulus(mag, distance)
    abs_mag_dist = distance_modulus(mag_dist, distance_dist)
    abs_mag_error= get_errors(abs_mag_dist)
    if plot_all:
        plt.hist(abs_mag_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
        plt.axvline(np.nanmedian(abs_mag_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
        plt.axvline(np.nanpercentile(abs_mag_dist, 84), color = 'cyan')
        plt.errorbar(abs_mag, 0.5, xerr = abs_mag_error, marker = '*', markersize = 8, color = 'b', label = "M_"+passband_string, capsize = 4)
        plt.xlabel('M_'+ passband_string)
        plt.legend()
        plt.show()
    else:
        pass
    return abs_mag, abs_mag_error, abs_mag_dist


def plot_target_table(input_table, absmag='g', colours= ['bp', 'rp'], list_color=list_color, pseudo_colour=False, annotate=annotate, label='', markersize=markersize,num="",use_primary_parallax=False):
    for row in input_table:
        if pseudo_colour:
            target_absmag, target_absmag_err, target_absmag_dist= get_pass_abs_mag(row, plot_all = False, passband_string= absmag)
            target_pseudo_colour=row['astrometric_pseudo_colour'+num]
            target_pseudo_dist= get_mc_distribution(target_pseudo_colour, row['astrometric_pseudo_colour_error'+num])
            target_pseudo_wave_dist=1./target_pseudo_dist*1e4
            target_colour_dif_err= get_errors(target_pseudo_wave_dist)
            target_colour_dif = 1./target_pseudo_colour*1e4
            plt.xlabel('1./astrometric_pseudo_colour (angstroms)')
            plt.errorbar(np.copy(target_colour_dif), np.copy(target_absmag), yerr = np.copy(target_absmag_err),  xerr = np.copy(target_colour_dif_err), marker = 'o', markersize = markersize, color = list_color, capsize = 4, linestyle ='none')
            print(target_colour_dif, target_absmag)
            print(row['name'+num])
            if annotate:
                plt.annotate(str(row['name'+num]),xy=(np.copy(target_colour_dif), np.copy(target_absmag)), xycoords='data', xytext=(np.copy(target_colour_dif+0.01),np.copy(target_absmag-0.1)), textcoords= 'data' , fontsize=8, color =list_color)
            else:
                pass
        else:
            print('===========')
            #print(row)
            target_absmag, target_absmag_err, target_absmag_dist= get_pass_abs_mag(row, plot_all = False, passband_string= absmag,num=num,use_primary_parallax=use_primary_parallax)
            target_colour_dif, target_colour_dif_err = get_colour_dif(row, plot_all = False, colours=colours,num=num)
            if row['phot_'+colours[0]+ '_mean_mag'+num]>1e18:
                target_colour_dif= x_fill
                target_colour_dif_err= 0
            if row['parallax'+num]>1e18:
                target_absmag=y_fill
            if error_bar:
                plt.errorbar(np.copy(target_colour_dif), np.copy(target_absmag), yerr = np.copy(target_absmag_err),  xerr = np.copy(target_colour_dif_err), marker = 'o', markersize = markersize, color = list_color, capsize = 4, linestyle ='none')
            else:
                plt.plot(np.copy(target_colour_dif), np.copy(target_absmag),  marker = 'o', markersize = markersize, color = list_color,  linestyle ='none')
            print(target_colour_dif, target_absmag)
            try:
                print(row['name'+num])
                if annotate:
                    plt.annotate(str(row['name'+num]),xy=(np.copy(target_colour_dif), np.copy(target_absmag)), xycoords='data', xytext=(np.copy(target_colour_dif+0.01),np.copy(target_absmag-0.1)), textcoords= 'data' , fontsize=8, color =list_color)
            except KeyError as newerror:
                print("KeyError:",newerror)
            else:
                pass
    if label=='':
        pass
    elif error_bar:
        plt.errorbar(np.nan, np.nan, yerr = np.nan,  xerr = np.nan, marker = 'o', markersize = markersize, color = list_color, capsize = 4, label =label, linestyle ='none')
    else:
        plt.plot(np.nan, np.nan,marker = 'o', markersize = markersize, color = list_color, label =label, linestyle ='none')
    return



def plot_bkg_cmd(generic_table= generic_table, absmag='g', colours=['bp','rp'], pseudo_colour= False):
    """
    Generates the hexbin histogram of whatever Gaia sample, such as the 100pc sample, or 47 Tucanae
    """
    
    try:
        if pseudo_colour:
            generic_parallax = generic_table ['parallax']+parallax_correction
            generic_parallax = generic_parallax *1e-3 #parallax in arcseconds now
            generic_distance = 1./generic_parallax #parsec distance
            generic_mag = generic_table['phot_'+absmag+'_mean_mag']
            generic_absmag = distance_modulus(generic_mag, generic_distance)
            generic_colour_dif = 1./generic_table['astrometric_pseudo_colour']*1e4
            plt.xlabel('1./astrometric_pseudo_colour (angstroms)')
        else:
            generic_parallax = generic_table ['parallax']+parallax_correction
            generic_parallax = generic_parallax *1e-3 #parallax in arcseconds now
            generic_distance = 1./generic_parallax #parsec distance
            generic_mag = generic_table['phot_'+absmag+'_mean_mag']
            generic_colour0 = generic_table['phot_'+colours[0]+'_mean_mag']
            generic_colour1= generic_table['phot_'+colours[1]+'_mean_mag']
            generic_colour_dif= generic_colour0- generic_colour1
            generic_absmag = distance_modulus(generic_mag, generic_distance)
            plt.xlabel(colours[0]+"-"+colours[1])
            plt.xlim(axes_x)


    except KeyError as error:
        print(error)
        print("assuming it's the simplified file.")
        if ((colours!=['bp','rp']) or (absmag != 'g')):
                print("Can't use the simplified file with alternative colours or absolute magnitudes; it only works with M_G and BP-RP")
        generic_absmag= generic_table['mg']
        generic_colour_dif= generic_table['bp_rp']
    #polything = plt.hexbin(generic_colour_dif, generic_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
    polything = plt.hexbin(generic_colour_dif, generic_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label='')
    counts = polything.get_array()
    print(counts.shape)
    counts= np.sqrt(counts)
    #counts=np.log(counts)
    polything.set_array(counts)
    polything.autoscale()
    plt.ylim(axes_y)
    plt.gca().invert_yaxis()
    #plt.xlabel(r'$G_{BP} - G_{RP}$')
    #plt.ylabel(r'$M_G$')
    plt.ylabel('M_'+absmag)
    plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)
    
    
    return

def fig_bkg_cmd(longax, generic_table= generic_table, absmag='g', colours=['bp','rp']):
    """
    Generates the hexbin histogram of whatever Gaia sample, such as the 100pc sample, or 47 Tucanae
    """
        
    try:
        generic_parallax = generic_table ['parallax']+parallax_correction
        generic_parallax = generic_parallax *1e-3 #parallax in arcseconds now
        generic_distance = 1./generic_parallax #parsec distance
        generic_mag = generic_table['phot_'+absmag+'_mean_mag']
        generic_colour0 = generic_table['phot_'+colours[0]+'_mean_mag']
        generic_colour1= generic_table['phot_'+colours[1]+'_mean_mag']
        generic_colour_dif= generic_colour0- generic_colour1
        generic_absmag = distance_modulus(generic_mag, generic_distance)

    except KeyError as error:
        print(error)
        print("assuming it's the simplified file.")
        if ((colours!=['bp','rp']) or (absmag != 'g')):
                print("Can't use the simplified file with alternative colours or absolute magnitudes; it only works with M_G and BP-RP")
        generic_absmag= generic_table['mg']
        generic_colour_dif= generic_table['bp_rp']
    polything = longax.hexbin(generic_colour_dif, generic_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1)
    counts = polything.get_array()
    print(counts.shape)
    #counts= np.sqrt(counts)
    counts=np.log(counts)
    polything.set_array(counts)
    polything.autoscale()
    longax.set_ylim(axes_y)
    plt.gca().invert_yaxis()
    #plt.xlabel(r'$G_{BP} - G_{RP}$')
    longax.set_xlabel(colours[0]+"-"+colours[1])
    longax.set_xlim(axes_x)
    #plt.ylabel(r'$M_G$')
    longax.set_ylabel('M_'+absmag)
    #plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)
    
    
    return longax


def make_cmd(target_table=target_table, generic_table= generic_table, absmag='g', colours=['bp', 'rp'], plot_cuts=False):
    plot_bkg_cmd(generic_table=generic_table, absmag=absmag, colours=colours)
    if plot_cuts:
        if ((colours==['bp','rp']) and (absmag=='g')):
            plot_nicola_cuts()
        elif((colours==['g','rp']) and (absmag=='g')):
            plot_ben_cuts()
        else:
            print('no cuts for selected absolute mag vs. colours')
    else:
        pass
    plot_target_table(target_table, absmag=absmag, colours=colours)
    plt.show()
    return


def plot_abs_v_abs(generic_table= generic_table, colours=['g','rp']):
    """
    generates the hexbin of two absolute magnitude plots
    """
    
    try:
        generic_parallax = generic_table ['parallax']+parallax_correction
        generic_parallax = generic_parallax *1e-3 #parallax in arcseconds now
        generic_distance = 1./generic_parallax #parsec distance
        mag0 = generic_table['phot_'+colours[0]+'_mean_mag']
        mag1=generic_table['phot_'+colours[1]+'_mean_mag']
        #generic_colour0 = generic_table['phot_'+colours[0]+'_mean_mag']
        #generic_colour1= generic_table['phot_'+colours[1]+'_mean_mag']
        #generic_colour_dif= generic_colour0- generic_colour1
        absmag0 = distance_modulus(mag0, generic_distance)
        absmag1=distance_modulus(mag1,generic_distance)
        plt.xlabel('M_'+colours[0])
        polything = plt.hexbin(absmag0, absmag1, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
        counts = polything.get_array()
        print(counts.shape)
        counts= np.sqrt(counts)
        #counts=np.log(counts)
        polything.set_array(counts)
        polything.autoscale()
        plt.ylim(axes_y)
        plt.gca().invert_yaxis()
        plt.xlim(axes_y)
        plt.gca().invert_xaxis()
        #plt.xlabel(r'$G_{BP} - G_{RP}$')
        #plt.ylabel(r'$M_G$')
        plt.ylabel('M_'+colours[1])
        plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)


    except KeyError as error:
        print(error)
        print("absolute mags generation failed")
    
    plt.show()
    
    return


#def make_cmd_from_list(target_list=[target_table], generic_table=generic_table, absmag='g', colours=['bp','rp']):
    #plot_bkg_cmd(generic_table=generic_table, absmag=absmag, colours=colours)
    #for target_table in target_list:
        #plot_target_table(target_table, absmag=absmag, colours=colours)
        

####################################

#longfig= plt.figure(figsize= (36, 36))
#longax= longfig.add_subplot(1,1,1)
#longax= fig_bkg_cmd(longax)
#longfig.tight_layout()
#longfig.savefig(str(distance)+'pc_big_cmd.png')
#longfig.closefig()
    #longax.plot(plot_waves, flux_normed, color= 'k')
    #plot_element_lines(plot_waves, longax)
    #plt.grid()
    ##longax.set_ylabel('Flux (normed)')
    #longax.set_ylabel('Flux (cgs units)')
    #longax.set_xlabel('Wavelength $(\AA)$')
    #longfig.savefig(dest_dir+target_dir+'long_spectrum_'+ str(wave_limits[0])+','+str(wave_limits[1]) + '.pdf')

if __name__ == '__main__':
    #plot_bkg_cmd()
    #plt.title(generic_input)
    #plt.show()
    #plot_bkg_cmd(absmag= absmag_band, colours= colours)
    #plt.title(generic_input)
    #plt.show()
    #plot_bkg_cmd(generic_table= generic_table, absmag= absmag_band, colours= ['bp','g'])
    #plt.title(generic_input)
    #plt.show()

    #plot_bkg_cmd(absmag= 'bp', colours= ['bp','g'])
    #plt.title(generic_input)
    #plt.show()


    #plot_bkg_cmd(absmag= 'rp', colours= ['g','rp'])
    #plt.title(generic_input)
    #plt.show()

    #plot_target_table(target_table, pseudo_colour=True)
    #plot_bkg_cmd(absmag='g', pseudo_colour=True)
    #plt.title(generic_input)
    #plt.show()
    #plot_abs_v_abs()
    #make_cmd(target_table=other_target_table, generic_table= generic_table)
    
    plot_bkg_cmd(generic_table=generic_table, colours=['g','rp'])
    plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=True)
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=annotate)
    #plot_target_table(target_table, colours=['g','rp'],num="",annotate=True)
    
    #down_mags=distance_modulus(target_table['phot_g_mean_mag'],1./(1e-3*target_table['parallax']))
    #down_inds=np.where(down_mags>12)
    #plot_target_table(target_table[down_inds], colours=['g','rp'],num="",annotate=True)
    plot_target_table(target_table, colours=['g','rp'],num="")
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='purple', annotate=False,label='No XP')
    #plt.title("WD+MS wide binary survey White Dwarfs")
    #plt.title('Vennes et al. 2023 new Li-polluted DZAH')
    #plt.title('WD J0212-5522 and its main-sequence companion')
    #plot_ben_cuts()
    plt.plot([1.55,-0.03],[16.75,4.58])
    plt.plot([0,1.3,1.3],[14.8,14.8,18])
    plt.legend()
    plt.show()
    
    plot_bkg_cmd(generic_table=generic_table, colours=['bp','rp'])
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=True)
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=annotate)
    #plot_target_table(target_table, colours=['g','rp'],num="2",list_color='g')
    plot_target_table(target_table, colours=['bp','rp'],num="",label='has XP')
    plot_target_table(other_target_table, colours=['bp','rp'], list_color='purple', annotate=False,label='No XP')
    plt.title("WD components of WD+MS wide binaries from Ben's 2023-09-05\ndown-select of El-Badry using the brightness-separation function")
    #plot_ben_cuts()
    #plt.plot([1.55,-0.03],[16.75,4.58])
    plt.legend()
    plt.show()
    
    plot_bkg_cmd(generic_table=generic_table, colours=['bp','rp'])
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=True)
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=annotate)
    #plot_target_table(target_table, colours=['g','rp'],num="2",list_color='g')
    plot_target_table(target_table, colours=['bp','rp'],num="")
    plot_target_table(other_target_table, colours=['bp','rp'], list_color='purple', annotate=True)
    #plot_ben_cuts()
    plt.plot([1.55,-0.03],[16.75,4.58])
    plt.show()


    plot_bkg_cmd(generic_table=generic_table, colours=['g','rp'])
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=True)
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=annotate)
    #plot_target_table(target_table, colours=['g','rp'],num="2",list_color='g')
    plot_target_table(target_table, colours=['g','rp'],num="2",list_color='g',use_primary_parallax=True)
    plot_target_table(target_table, colours=['g','rp'],num="1")
    plot_target_table(other_target_table, colours=['g','rp'], list_color='purple', annotate=True)
    #plot_ben_cuts()
    plt.plot([1.55,-0.03],[16.75,4.58])
    plt.show()
    
    plot_bkg_cmd(generic_table=generic_table, colours=['bp','rp'])
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=True)
    #plot_target_table(other_target_table, colours=['g','rp'], list_color='r', annotate=annotate)
    plot_target_table(target_table, colours=['bp','rp'],num="1")
    plot_target_table(target_table, colours=['bp','rp'],num="2",list_color='g')
    plt.show()


    plot_bkg_cmd(generic_table=generic_table, colours=['g','rp'], absmag='rp')
    plot_target_table(other_target_table, colours=['g','rp'], absmag='rp', list_color='r', annotate=True)
    plot_target_table(target_table, colours=['g','rp'], absmag='rp')
    plt.show()


    plot_bkg_cmd(generic_table=generic_table, colours=['g','rp'], absmag='bp')
    plot_target_table(other_target_table, colours=['g','rp'], absmag='bp', list_color='r', annotate=True)
    plot_target_table(target_table, colours=['g','rp'], absmag='bp')
    plt.show()


    plot_target_table(other_target_table, colours=['bp','rp'], list_color='r', annotate=True)
    plot_target_table(target_table, colours=['bp','rp'])
    plot_bkg_cmd(generic_table=generic_table, colours=['bp','rp'])
    plt.show()


    #plot_target_table(target_table, colours=['bp','g'])
    #plot_target_table(other_target_table, colours=['bp','g'], list_color='r', annotate=True)
    #plot_bkg_cmd(generic_table=generic_table, colours=['bp','g'])
    #plt.show()


    #plot_bkg_cmd(generic_table=generic_table, absmag= 'rp', colours=['bp','rp'])
    #plot_target_table(target_table, colours=['bp','rp'], absmag= 'rp')
    #plot_target_table(other_target_table, absmag= 'rp', colours=['bp','rp'], list_color='r')
    #plt.show()


    #plot_bkg_cmd(generic_table=generic_table, absmag= 'bp', colours=['bp','rp'])
    #plot_target_table(target_table, colours=['bp','rp'], absmag= 'bp')
    #plot_target_table(other_target_table, absmag= 'bp', colours=['bp','rp'], list_color='r')
    #plt.show()

    #plot_bkg_cmd(generic_table=generic_table, colours=['bp','rp'])
    #plot_target_table(target_table, colours=['bp','rp'])
    #plot_target_table(other_target_table, colours=['bp','rp'], list_color='r')
    #plt.show()

    #make_cmd(target_table=target_table, generic_table= generic_table)
    make_cmd(target_table=target_table, generic_table= generic_table, plot_cuts=True)
    make_cmd(target_table=target_table, generic_table= generic_table, absmag= absmag_band, colours= colours, plot_cuts=True)
    make_cmd(target_table=target_table, generic_table= generic_table, absmag= absmag_band, colours= ['bp','g'])
    #make_cmd(target_table=target_table, generic_table= generic_table, absmag= 'bp', colours= ['bp','rp'])

    #make_cmd(target_table=target_table, generic_table= generic_table, absmag= 'bp', colours= ['bp','rp'])
