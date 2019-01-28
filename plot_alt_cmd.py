"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-01-07

This script should plot lists of targets over the CMD for a given distance, and allow for the creation of different
CMDs than M_G vs. BP-RP, i.e. M_G vs. G-RP

It's going to borrow heavily from plot_gaia_cmd.py, but I'm going to try to keep it clean of all of the extra stuff
I have in there to do calculations of a single target.

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





axes_x= [-1.5, 6]
axes_y = [-4,18]


#axes_x= [-3, 6]
#axes_y = [-5,25]

y_fill= axes_y[1]
x_fill= axes_x[0]

list_color = '#1ca1f2'
single_list=False#turns off the original for-loop method for plotting from a single list
error_bar=True #turns off error bars on the multiple list plot, meaning it has no effect on anything if single_list==True
annotate= False #controls whether or not object names appear beside points in the scatter plots. Should be turned off for >~20 targets appearing close together
parallax_correction = 0.029 #from Lindgren et al 2018


#######3error distribution variables
mc_number = 10000
percent_off = 34 #1-sigma equivalent
#############

target_input='20190107_chris_merge_gaia.csv'
target_input=  'Lindegren_appC_selB_antiC_cut2_gaia_sc.csv'
#target_input='l-0.3bp_g_gaia_corr_full.csv'
#target_input='elm_survey_gaia.csv'
#target_input='pre_elms_gaia.csv'
#target_input = '20190109_blue_gaia.csv'
#target_input='hot_wind_wds_gaia.csv'
#target_input='Eriks_disk_candidates_gaia.csv'
#target_input='dC_sample_roulston2018_gaia.csv'
#target_input='20190111_red_things_gaia.csv'
#target_input='alt_red_things_g_rp_greater_17_gaia.csv'
#target_input='WD_cooling_tip_gaia.csv'
#target_input='weird_CPM_binary_gaia.csv'
#target_input='20190121_excess_interesting_gaia.csv'
#target_input= 'ar_sco_gaia.csv'
#target_input='20190123_new_red_things_gaia_sc.csv'
#target_input= 'DQpec_gaia.csv'
target_label= ''

#num_targs = 'all'
#num_targs = '47Tuc'
num_targs= 'Lindegren'
selection_letter= 'C'
distance = 200
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

####################################
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
    except astropy.units.core.UnitsError as error:
        return  np.array([[median.value-low_bar],[high_bar-median.value]])


def get_filter_vals(table, filter_string):
    flux_string = 'phot_'+filter_string+ '_mean_flux'
    phot_mean_flux = table[flux_string]
    error_string = flux_string + '_error'
    phot_mean_flux_error = table[error_string]
    flux_distribution = get_mc_distribution(phot_mean_flux, phot_mean_flux_error)
    return phot_mean_flux, flux_distribution



def get_colour_dif(table, plot_all = False, verbose =True, colours=['bp','rp']):
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
    mean_flux0, dist0 = get_filter_vals(table, colours[0])
    mean_flux1  , dist1 = get_filter_vals(table, colours[1])
    mag0 = get_mag(mean_flux0, colours[0])
    mag1 = get_mag(mean_flux1, colours[1])
    if verbose:
        print(colours[0]+"_calc-bp_measured", mag0 - table['phot_' + colours[0]+ '_mean_mag'])
        print(colours[1]+ "_calc - rp_measured", mag1 - table['phot_' + colours[1]+'_mean_mag'])
    mag_dist0 = get_mag(dist0, colours[0])
    mag_dist1 = get_mag(dist1, colours[1])
    colour_dif = mag0- mag1
    colour_dif_dist= mag_dist0- mag_dist1
    colour_dif_error = get_errors(colour_dif_dist)
    if plot_all:
        plt.hist(bp_rp_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
        plt.axvline(np.nanmedian(colour_dif_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
        plt.axvline(np.nanpercentile(colour_dif_dist, 84), color = 'cyan')
        plt.errorbar(colour_dif, 0.5, xerr = colour_dif_error, marker = '*', markersize = 8, color = 'b', label = colours[0]+"-" + colours[1], capsize = 4)
        #plt.xlabel(r'$G_{BP}-G_{RP}$')
        plt.xlabel(colours[0]+'-'+colours[1])
        plt.legend()
        plt.show()
    return colour_dif, colour_dif_error




def get_pass_abs_mag(table, plot_all = False, passband_string= 'g', verbose = True):
    mean_flux, flux_dist = get_filter_vals(table, passband_string)
    mag = get_mag(mean_flux, passband_string)
    flux_dist= remove_negative(flux_dist, verbose=verbose)
    mag_dist= get_mag(flux_dist, passband_string)
    parallax = table['parallax']+parallax_correction
    parallax = parallax*1e-3
    distance = 1./parallax
    parallax_error = table['parallax_error']*1e-3
    parallax_dist = get_mc_distribution(parallax, parallax_error)
    parallax_dist = remove_negative(parallax_dist, verbose= verbose)
    if verbose:
        print(passband_string+ "_calc" + "-" + passband_string+ "_measured", mag - table['phot_' +passband_string+'_mean_mag'])
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


def plot_target_table(input_table, absmag='g', colours= ['bp', 'rp']):
    for row in input_table:
        print('===========')
        #print(row)
        target_absmag, target_absmag_err, target_absmag_dist= get_pass_abs_mag(row, plot_all = False, passband_string= absmag)
        target_colour_dif, target_colour_dif_err = get_colour_dif(row, plot_all = False, colours=colours)
        if row['phot_'+colours[0]+ '_mean_mag']>1e18:
            target_colour_dif= x_fill
            target_colour_dif_err= 0
        if row['parallax']>1e18:
            target_absmag=y_fill
        plt.errorbar(np.copy(target_colour_dif), np.copy(target_absmag), yerr = np.copy(target_absmag_err),  xerr = np.copy(target_colour_dif_err), marker = 'o', markersize = 6, color = list_color, capsize = 4, label = target_label, linestyle ='none')
        print(target_colour_dif, target_absmag)
        print(row['name'])
        if annotate:
            plt.annotate(str(row['name']),xy=(np.copy(target_colour_dif), np.copy(target_absmag)), xycoords='data', xytext=(np.copy(target_colour_dif+0.01),np.copy(target_absmag-0.1)), textcoords= 'data' , fontsize=8, color =list_color)
        else:
            pass
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
    polything = plt.hexbin(generic_colour_dif, generic_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
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
    polything = longax.hexbin(generic_colour_dif, generic_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
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


def make_cmd(target_table=target_table, generic_table= generic_table, absmag='g', colours=['bp', 'rp']):
    plot_bkg_cmd(generic_table=generic_table, absmag=absmag, colours=colours)
    plot_target_table(target_table, absmag=absmag, colours=colours)
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

#plot_bkg_cmd(absmag='g', pseudo_colour=True)
#plt.title(generic_input)
#plt.show()

make_cmd(target_table=target_table, generic_table= generic_table)
make_cmd(target_table=target_table, generic_table= generic_table, absmag= absmag_band, colours= colours)
make_cmd(target_table=target_table, generic_table= generic_table, absmag= absmag_band, colours= ['bp','g'])
