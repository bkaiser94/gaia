"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-01-09

This should combine star magnitudes in the Gaia CMD to simulate an unresolved binary's position in the diagram with an empirical method.



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

import gaia_extinction
#import wdatmos
import plotting_dicts as pod

#merged_file='WISEA0615m1247.csv'

#star1_file='WDJ2356m209.csv'

#star 2 is the one that should be solved for using the file.


absmag_band= 'g'
colours= ['g','rp']





axes_x= [-1.5, 6]
axes_y = [-4,18]

y_fill= axes_y[1]
x_fill= axes_x[0]

star1_color = '#1ca1f2'
list_color=star1_color
merged_color= 'magenta'
star2_color= 'r'
single_list=True#turns off the original for-loop method for plotting from a single list
error_bar=True #turns off error bars on the multiple list plot, meaning it has no effect on anything if single_list==True
annotate= True #controls whether or not object names appear beside points in the scatter plots. Should be turned off for >~20 targets appearing close together
parallax_correction = 0.029 #from Lindgren et al 2018


#######3error distribution variables
mc_number = 10000
percent_off = 34 #1-sigma equivalent
#############
num_targs='all'
#num_targs = 'Lindegren'
#num_targs = '47Tuc'
distance = 200
grid_num = 225
selection_letter= 'B'

#num_stars = 20 #number of stars in the track to use for merging
num_stars = 5 #number of stars in the track to use for merging

extrap_dist = 0.0 #magnitude distance to go outside the line ends for the hypothetical star mergers
extrap_dist1= 0.1
extrap_dist2=extrap_dist
#############3

#star1_input= '20190109_star1_ultracool.csv'
#merged_star_input= '20190109_merged_star.csv'


#star1_input= '20190109_star1_Erik.csv'
#merged_star_input= '20190109_merged_star_Erik.csv'

#star1_input= 'weird_CPM_binary_gaia.csv'

#star1_input= 'test_WD.csv'
#star1_input= 'mystery_red_object.csv'
star1_input='2MASSJ1458p2839_gaia.csv'
#star1_input='2MASSJ1055p0808_gaia.csv'
#star1_input= '2MASSIJ0821p4532_gaia.csv'
#star1_input='sdssj1408p2021_gaia.csv'

#star1_input='two_low_reds_gaia.csv'
#star2_input= '20190218_test_ultcool.csv'
star2_input='sdssj1330p6435_gaia.csv'
#star1_input= 'CE40_gaia.csv'
#merged_star_input= 'WISEA0615m1247.csv'

merged_star_input='WISEA0615m1247.csv'
#merged_star_input='LP178m49_gaia.csv' 
#merged_star_input= 'WD1133p358_gaia.csv'
#merged_star_input='WISEA0238p3617.csv'

#star1_input='sdssj1246p3608_gaia.csv'
#star1_input='sdssj1408p2021_gaia.csv'


target_label = ''


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
    #generic_input='Lindegren_appC_selC_noBLMC.csv'
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
#target_table = Table.read(target_input)
star1_table= Table.read(star1_input)[0]
merged_star_table= Table.read(merged_star_input)[0]
star2_table= Table.read(star2_input)[0]
input_table= Table.read(star1_input)
print(star1_table)
print(merged_star_table)

####################################
def distance_modulus(g_mag, distance):
    return g_mag - 5*np.log10(distance/10.)
    
def get_mag(flux, filter_string):
    mag0 = zeropoint_dict[filter_string][0]
    return -2.5*np.log10(flux) +mag0

def mag_to_flux(mag, filter_string):
    mag0=zeropoint_dict[filter_string][0]
    return 10.**((mag-mag0)/-2.5)

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
        print(change_array[:min_inds], match_array[:min_inds])
        return change_array[:min_inds], match_array[:min_inds]
    except AttributeError:
        #the inputs aren't actually arrays
        print("inputs aren't actually arrays")
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

def get_star_absmags(star_table, plot_all=False):
    g_absmag, g_absmag_error, g_absmag_dist= get_pass_abs_mag(star_table, passband_string='g', plot_all= plot_all)
    bp_absmag, bp_absmag_error, bp_absmag_dist= get_pass_abs_mag(star_table, passband_string='bp', plot_all=plot_all)
    rp_absmag, rp_absmag_error, rp_absmag_dist= get_pass_abs_mag(star_table, passband_string='rp', plot_all=plot_all)
    star_dict= {'g':{
        'absmag':g_absmag,
        'absmag_error':g_absmag_error,
        'absmag_dist':g_absmag_dist},
    'bp':{
        'absmag':bp_absmag,
        'absmag_error':bp_absmag_error,
        'absmag_dist':bp_absmag_dist},
    'rp':{
        'absmag':rp_absmag,
        'absmag_error':rp_absmag_error,
        'absmag_dist':rp_absmag_dist}}
    return star_dict

def generate_star_track(endpoint1,endpoint2,mode='linear', num_points=num_stars, colours= ['g','rp'], extrap_dist= 0.0):
    """
    Takes two CMD endpoints and draws a line between them then populates with evenly-spaced stars (in
    magnitude space).
    
    The inputs will be bp-rp I think... I guess I'll need to make this flexible in the near future.
    Currently it actually has to be G-RP because that's the color that can be affixed for combination
    with M_G affixed.
    
    I'll output a list of dicts I guess that can be accessed for each star
    """
    xpoints= np.array([endpoint1[0], endpoint2[0]])
    ypoints= np.array([endpoint1[1],endpoint2[1]])
    if mode=='linear':
        coeffs = np.polyfit(xpoints, ypoints, 1)
    else:
        pass
    
    star_track_x = np.linspace(xpoints[0]-extrap_dist,xpoints[1]+extrap_dist,num_points)
    star_track_y=np.polyval(coeffs, star_track_x)
    star_track_g_absmag= star_track_y
    star_track_rp_absmag= star_track_y-star_track_x
    dict_list= []
    for index in range(0,num_points):
        star_dict= {'g':{
            'absmag':star_track_g_absmag[index],
            'absmag_error': 0,
            'absmag_dist':0},
        'bp':{
            'absmag':0,
            'absmag_error':0,
            'absmag_dist':0},
        'rp':{
            'absmag':star_track_rp_absmag[index],
            'absmag_error':0,
            'absmag_dist':0}}
        dict_list.append(star_dict)
    return dict_list

def generate_star_grid(endpoint1,endpoint2,gridrows= 4, grid_space= 0.05, mode='linear', num_points=num_stars, colours= ['g','rp'], extrap_dist= 0.0):
    
    return

def generate_random_grid(abs_bounds=[14.,17.], colour_bounds=[-0.5, 1.], num_points=num_stars):
    """
    Create a grid of stars randomly distributed in the region defined to be used either as merged stars  or 
    components in the other functions.
    
    """
    abs_mags = np.random.rand(num_points)*(abs_bounds[1]-abs_bounds[0])+abs_bounds[0]
    colour_vals= np.random.rand(num_points)*(colour_bounds[1]-colour_bounds[0])+colour_bounds[0]
    rp_vals= abs_mags-colour_vals
    print("abs_mags", abs_mags)
    print("colour_vals", colour_vals)
    print("rp_vals", rp_vals)
    dict_list= []
    for index in range(0,num_points):
        star_dict= {'g':{
            'absmag':abs_mags[index],
            'absmag_error': 0,
            'absmag_dist':0},
        'bp':{
            'absmag':0,
            'absmag_error':0,
            'absmag_dist':0},
        'rp':{
            'absmag':rp_vals[index],
            'absmag_error':0,
            'absmag_dist':0}}
        dict_list.append(star_dict)
    for thing in dict_list:
        print(thing)
    return dict_list


def find_star2(star1_table, merged_star_table, bounded=False, bounds=[], halfreal=False, plot_all=False, num_points=num_stars):
    """
    Take the table for the first star, and the star that represents the position of the end product, and then
    return the position of the second star that is required to end up in that spot.
    
    """
    if halfreal:
        g_rp_vals= []
        m_g_vals=[]
        merged_absmag_dict =get_star_absmags(merged_star_table, plot_all=False)
        #star_track1= generate_star_track([1.22, 11.],[1.56,15.226], num_points= num_points, extrap_dist=extrap_dist1)
        star_track1=generate_random_grid(abs_bounds=[merged_absmag_dict['g']['absmag'], 17.])
        merged_g_flux= mag_to_flux(merged_absmag_dict['g']['absmag'],'g')
        merged_rp_flux= mag_to_flux(merged_absmag_dict['rp']['absmag'],'rp')
        g_rp_subvals=[]
        m_g_subvals=[]
        plt.errorbar(merged_absmag_dict['g']['absmag']- merged_absmag_dict['rp']['absmag'], merged_absmag_dict['g']['absmag'], yerr = merged_absmag_dict['g']['absmag_error'],  xerr = get_errors(merged_absmag_dict['g']['absmag_dist']-merged_absmag_dict['rp']['absmag_dist']), marker = 'o', markersize = 6, color =merged_color, capsize = 4, label = 'merged', linestyle ='none')
        star_track1_mg=[]
        star_track1_g_rp=[]
        for star1 in star_track1:
            star1_g_flux= mag_to_flux(star1['g']['absmag'],'g')
            star1_rp_flux= mag_to_flux(star1['rp']['absmag'],'rp')
            #plt.plot(star1['g']['absmag']-star1['rp']['absmag'], star1['g']['absmag'], color= 'g', linestyle= None, marker='o', markersize= 6)
            star_track1_mg.append(star1['g']['absmag'])
            star_track1_g_rp.append(star1['g']['absmag']-star1['rp']['absmag'])
            star2_g_flux = merged_g_flux-star1_g_flux
            star2_rp_flux= merged_rp_flux-star1_rp_flux
            star2_g_absmag= get_mag(star2_g_flux, 'g')
            star2_rp_absmag= get_mag(star2_rp_flux, 'rp')
            g_rp_subvals.append(star2_g_absmag-star2_rp_absmag)
            m_g_subvals.append(star2_g_absmag)
            g_rp_vals.append(star2_g_absmag-star2_rp_absmag)
            m_g_vals.append(star2_g_absmag)
            #plt.plot(g_rp_subvals, m_g_subvals, color='magenta', linestyle= None, marker='o', markersize= 6)
        number_list=[]
        plt.plot(g_rp_subvals, m_g_subvals, color=star2_color, linestyle= None, marker='o', markersize= 6)
        plt.plot(star_track1_g_rp, star_track1_mg, color=star1_color, linestyle= None, marker='o', markersize= 6)
        
        #plt.annotate(str(),xy=(np.copy(target_colour_dif), np.copy(target_absmag)), xycoords='data', xytext=(np.copy(target_colour_dif+0.01),np.copy(target_absmag-0.1)), textcoords= 'data' , fontsize=8, color =list_color)
        return g_rp_vals, m_g_vals
    else:
        merged_absmag_dict= get_star_absmags(merged_star_table, plot_all=False)
        
        def get_star2_absmag(bandpass='g'):
            try:
                print(bandpass, 'mag:', merged_absmag_dict[bandpass]['absmag'])
                star2_flux= mag_to_flux(merged_absmag_dict[bandpass]['absmag'],bandpass)-mag_to_flux(star1_absmag_dict[bandpass]['absmag'], bandpass)
                print(bandpass, 'flux:', star2_flux)
                star2_flux_dist= mag_to_flux(merged_absmag_dict[bandpass]['absmag_dist'],bandpass)-mag_to_flux(star1_absmag_dict[bandpass]['absmag_dist'], bandpass)
                star2_absmag= get_mag(star2_flux,bandpass)
                star2_absmag_dist= get_mag(star2_flux_dist,bandpass)
            except KeyError:
                print(bandpass, 'mag:', merged_absmag_dict[bandpass]['absmag'])
                star2_flux= mag_to_flux(merged_absmag_dict[bandpass]['absmag'],bandpass)-mag_to_flux(star1_absmag_dict[bandpass]['absmag'], bandpass)
                print(bandpass, 'flux:', star2_flux)
                star2_flux_dist= mag_to_flux(merged_absmag_dict[bandpass]['absmag_dist'],bandpass)-mag_to_flux(star1_absmag_dict[bandpass]['absmag'], bandpass)
                star2_absmag= get_mag(star2_flux,bandpass)
                star2_absmag_dist= get_mag(star2_flux_dist,bandpass)
            return star2_absmag, star2_absmag_dist
        star2_g_absmag, star2_g_absmag_dist= get_star2_absmag(bandpass='g')
        star2_bp_absmag, star2_bp_absmag_dist= get_star2_absmag(bandpass='bp')
        star2_rp_absmag, star2_rp_absmag_dist= get_star2_absmag(bandpass='rp')
        star2_bp_absmag_dist=remove_negative(star2_bp_absmag_dist)
        star2_rp_absmag_dist= remove_negative(star2_rp_absmag_dist)
        star2_rp_absmag_dist,star2_bp_absmag_dist= match_sizes(star2_rp_absmag_dist, star2_bp_absmag_dist)
        star2_g_absmag_dist,star2_rp_absmag_dist= match_sizes(star2_g_absmag_dist, star2_rp_absmag_dist)
        if plot_all:
            plt.title('')
            plt.show()
            
            print(star2_rp_absmag_dist)
            plt.title('star2 rp')
            plt.hist(star2_rp_absmag_dist)
            plt.show()
            
            plt.title('star 2 g')
            plt.hist(star2_g_absmag_dist)
            plt.show()
            
            plt.title('star2 bp')
            plt.hist(star2_bp_absmag_dist)
            plt.show()
        else:
            pass
        
        
        star2_bp_rp= star2_bp_absmag-star2_rp_absmag
        star2_bp_rp_dist= star2_bp_absmag_dist-star2_rp_absmag_dist
        star2_bp_rp_error=get_errors(star2_bp_rp_dist)
        star2_g_absmag_error= get_errors(star2_g_absmag_dist)
        
        star2_g_rp= star2_g_absmag-star2_rp_absmag
        star2_g_rp_dist= star2_g_absmag_dist-star2_rp_absmag_dist
        star2_g_rp_error=get_errors(star2_g_rp_dist)
        #star2_g_absmag_error= get_errors(star2_g_absmag_dist)
    
        #plt.title('bp_rp')
        #plt.hist(star2_bp_rp_dist)
        #plt.show()
        
        star1_bp_rp, star1_bp_rp_error = get_colour_dif(star1_table, colours=['bp','rp'])
        merged_bp_rp, merged_bp_rp_error= get_colour_dif(merged_star_table, colours=['bp','rp'])
        
        star1_g_rp, star1_g_rp_error = get_colour_dif(star1_table, colours=['g','rp'])
        merged_g_rp, merged_g_rp_error= get_colour_dif(merged_star_table, colours=['g','rp'])
        print('=======')
        print('star2 M_G:', star2_g_absmag, '+/-', star2_g_absmag_error)
        print('star2 BP-RP:', star2_bp_rp, '+/-', star2_bp_rp_error)
        #plt.errorbar(star2_bp_rp, star2_g_absmag, yerr = star2_g_absmag_error,  xerr = star2_bp_rp_error, marker = 'o', markersize = 6, color = 'magenta', capsize = 4, label = 'star2', linestyle ='none')
        #plt.errorbar(star1_bp_rp, star1_absmag_dict['g']['absmag'], yerr = star1_absmag_dict['g']['absmag_error'],  xerr = star1_bp_rp_error, marker = 'o', markersize = 6, color = list_color, capsize = 4, label = 'star1', linestyle ='none')
        #plt.errorbar(merged_bp_rp, merged_absmag_dict['g']['absmag'], yerr = merged_absmag_dict['g']['absmag_error'],  xerr = merged_bp_rp_error, marker = 'o', markersize = 6, color = 'g', capsize = 4, label = 'merged_star', linestyle ='none')
        plt.errorbar(star2_g_rp, star2_g_absmag, yerr = star2_g_absmag_error,  xerr = star2_g_rp_error, marker = 'o', markersize = 6, color = 'magenta', capsize = 4, label = 'star2', linestyle ='none')
        plt.errorbar(star1_g_rp, star1_absmag_dict['g']['absmag'], yerr = star1_absmag_dict['g']['absmag_error'],  xerr = star1_g_rp_error, marker = 'o', markersize = 6, color = list_color, capsize = 4, label = 'star1', linestyle ='none')
        plt.errorbar(merged_g_rp, merged_absmag_dict['g']['absmag'], yerr = merged_absmag_dict['g']['absmag_error'],  xerr = merged_g_rp_error, marker = 'o', markersize = 6, color = 'g', capsize = 4, label = 'merged_star', linestyle ='none')
        star2_dict= {'g':{
            'absmag':star2_g_absmag,
            'absmag_error':star2_g_absmag_error,
            'absmag_dist':star2_g_absmag_dist},
        'bp':{
            'absmag':star2_bp_absmag,
            'absmag_error':get_errors(star2_bp_absmag_dist),
            'absmag_dist':star2_bp_absmag_dist},
        'rp':{
            'absmag':star2_rp_absmag,
            'absmag_error':get_errors(star2_rp_absmag_dist),
            'absmag_dist':star2_rp_absmag_dist}}
    
        return star2_dict


def merge_stars(star1_table=[], star2_table=[], real_stars=False, bounded=False, bounds=[], num_points=num_stars, halfreal=False):
    """
   Take 2 stars and combine their fluxes to plot a new star on the CMD. This should be flexible enough to work
   with cooling sequences... although I guess I could just loop through it a bunch of times. Yeah, I guess
   I should just loop it externally.
    
    """
    #this is for the merging of two actual stars. I should probably have been calling this blending, in hindsight
    #this will need to be a try-except statement so it works with non-Gaia data inputs.
    if real_stars:
        star1_absmag_dict =get_star_absmags(star1_table, plot_all=True)
        star2_absmag_dict= get_star_absmags(star2_table, plot_all=True)
        def get_merged_absmag(bandpass='g'):
            print(bandpass, 'mag:', star2_absmag_dict[bandpass]['absmag'])
            merged_flux= mag_to_flux(star2_absmag_dict[bandpass]['absmag'],bandpass)+mag_to_flux(star1_absmag_dict[bandpass]['absmag'], bandpass)
            print(bandpass, 'flux:', merged_flux)
            merged_flux_dist= mag_to_flux(star2_absmag_dict[bandpass]['absmag_dist'],bandpass)+mag_to_flux(star1_absmag_dict[bandpass]['absmag_dist'], bandpass)
            merged_absmag= get_mag(merged_flux,bandpass)
            merged_absmag_dist= get_mag(merged_flux_dist,bandpass)
            return merged_absmag, merged_absmag_dist
        merged_g_absmag, merged_g_absmag_dist= get_merged_absmag(bandpass='g')
        merged_bp_absmag, merged_bp_absmag_dist= get_merged_absmag(bandpass='bp')
        merged_rp_absmag, merged_rp_absmag_dist= get_merged_absmag(bandpass='rp')
        merged_bp_absmag_dist=remove_negative(merged_bp_absmag_dist)
        merged_rp_absmag_dist= remove_negative(merged_rp_absmag_dist)
        merged_rp_absmag_dist,merged_bp_absmag_dist= match_sizes(merged_rp_absmag_dist, merged_bp_absmag_dist)
        
        merged_rp_absmag_dist,merged_g_absmag_dist= match_sizes(merged_rp_absmag_dist, merged_g_absmag_dist)
        
        plt.title('')
        plt.legend()
        plt.show()
        
        print(merged_rp_absmag_dist)
        plt.title('rp')
        plt.hist(merged_rp_absmag_dist)
        plt.show()
        
        plt.title('g')
        plt.hist(merged_g_absmag_dist)
        plt.show()
        
        plt.title('bp')
        plt.hist(merged_bp_absmag_dist)
        plt.show()
        
        
        #merged_bp_rp= merged_bp_absmag-merged_rp_absmag
        #merged_bp_rp_dist= merged_bp_absmag_dist-merged_rp_absmag_dist
        #merged_bp_rp_error=get_errors(merged_bp_rp_dist)
        #merged_g_absmag_error= get_errors(merged_g_absmag_dist)
        
        merged_g_rp= merged_g_absmag-merged_rp_absmag
        merged_g_rp_dist= merged_g_absmag_dist-merged_rp_absmag_dist
        merged_g_rp_error=get_errors(merged_g_rp_dist)
        merged_g_absmag_error= get_errors(merged_g_absmag_dist)
    
        #plt.title('bp_rp')
        #plt.hist(merged_bp_rp_dist)
        #plt.show()
        
        star1_g_rp, star1_g_rp_error = get_colour_dif(star1_table, colours=['g','rp'])
        star2_g_rp, star2_g_rp_error= get_colour_dif(star2_table, colours=['g','rp'])
        #merged_bp_rp, merged_bp_rp_error= get_colour_dif(merged_star_table, colours=['bp','rp'])
        
        print('=======')
        print('merged M_G:', merged_g_absmag, '+/-', merged_g_absmag_error)
        print('merged G-RP:', merged_g_rp, '+/-', merged_g_rp_error)
        plt.errorbar(merged_g_rp, merged_g_absmag, yerr = merged_g_absmag_error,  xerr = merged_g_rp_error, marker = 'o', markersize = 6, color = 'magenta', capsize = 4, label = 'merged', linestyle ='none')
        plt.errorbar(star1_g_rp, star1_absmag_dict['g']['absmag'], yerr = star1_absmag_dict['g']['absmag_error'],  xerr = star1_g_rp_error, marker = 'o', markersize = 6, color = list_color, capsize = 4, label = 'star1', linestyle ='none')
        plt.errorbar(star2_g_rp, star2_absmag_dict['g']['absmag'], yerr = star2_absmag_dict['g']['absmag_error'],  xerr = star2_g_rp_error, marker = 'o', markersize = 6, color = 'g', capsize = 4, label = 'star2', linestyle ='none')
        #pass
    elif halfreal:
        g_rp_vals= []
        m_g_vals=[]
        star1_absmag_dict =get_star_absmags(star1_table, plot_all=True)
        star_track2= generate_star_track([1.22, 11.],[1.56,15.226], num_points= num_points, extrap_dist=extrap_dist1)
        
        star1_g_flux= mag_to_flux(star1_absmag_dict['g']['absmag'],'g')
        star1_rp_flux= mag_to_flux(star1_absmag_dict['rp']['absmag'],'rp')
        g_rp_subvals=[]
        m_g_subvals=[]
        plt.errorbar(star1_absmag_dict['g']['absmag']-star1_absmag_dict['rp']['absmag'], star1_absmag_dict['g']['absmag'], yerr = star1_absmag_dict['g']['absmag_error'],  xerr = get_errors(star1_absmag_dict['g']['absmag_dist']-star1_absmag_dict['rp']['absmag_dist']), marker = 'o', markersize = 6, color = list_color, capsize = 4, label = 'star1', linestyle ='none')
        for star2 in star_track2:
            star2_g_flux= mag_to_flux(star2['g']['absmag'],'g')
            star2_rp_flux= mag_to_flux(star2['rp']['absmag'],'rp')
            plt.plot(star2['g']['absmag']-star2['rp']['absmag'], star2['g']['absmag'], color= 'g', linestyle= None, marker='o', markersize= 6)
            merged_g_flux = star1_g_flux+star2_g_flux
            merged_rp_flux= star1_rp_flux+star2_rp_flux
            merged_g_absmag= get_mag(merged_g_flux, 'g')
            merged_rp_absmag= get_mag(merged_rp_flux, 'rp')
            g_rp_subvals.append(merged_g_absmag-merged_rp_absmag)
            m_g_subvals.append(merged_g_absmag)
            g_rp_vals.append(merged_g_absmag-merged_rp_absmag)
            m_g_vals.append(merged_g_absmag)
            plt.plot(g_rp_subvals, m_g_subvals, color='magenta', linestyle= None, marker='o', markersize= 6)
        return g_rp_vals, m_g_vals
    else:
        #for star tracks
        print("has to be G-RP colors")
        
        #star_track2= generate_star_track([-0.127, 11.05], [0.8,15.5],num_points=num_points, extrap_dist=0.1) #approximate WD track
        #star_track2= generate_star_track([-0.16, 16.44], [0.84,15.61],num_points=num_points, extrap_dist=0.1) #approximate ultracool WD track
        star_track2= generate_star_track([1.22, 11.],[1.56,15.226], num_points= num_points, extrap_dist=extrap_dist1) #approximate Mstars
        star_track1= generate_star_track([1.22, 11.],[1.56,15.226], num_points= num_points, extrap_dist=extrap_dist1) #approximate Mstars
        #star_track1= generate_star_track([0.524, 5.178],[0.949,8.104], num_points= num_points, extrap_dist=extrap_dist2) #approximate main-sequence A-G
        #star_track2= generate_star_track([0.527, 5.28],[0.89,7.75], num_points= num_points) #approximate Mstars

        g_rp_vals= []
        m_g_vals=[]
        for star1 in star_track1:
            plt.plot(star1['g']['absmag']-star1['rp']['absmag'], star1['g']['absmag'],color='r', linestyle=None, marker='o', markersize= 6)
            star1_g_flux= mag_to_flux(star1['g']['absmag'],'g')
            star1_rp_flux= mag_to_flux(star1['rp']['absmag'],'rp')
            g_rp_subvals=[]
            m_g_subvals=[]
            for star2 in star_track2:
                plt.plot(star2['g']['absmag']-star2['rp']['absmag'], star2['g']['absmag'], color= '#1ca1f2', linestyle= None, marker='o', markersize= 6)
                star2_g_flux= mag_to_flux(star2['g']['absmag'],'g')
                star2_rp_flux= mag_to_flux(star2['rp']['absmag'],'rp')
                merged_g_flux = star1_g_flux+star2_g_flux
                merged_rp_flux= star1_rp_flux+star2_rp_flux
                merged_g_absmag= get_mag(merged_g_flux, 'g')
                merged_rp_absmag= get_mag(merged_rp_flux, 'rp')
                g_rp_vals.append(merged_g_absmag-merged_rp_absmag)
                m_g_vals.append(merged_g_absmag)
                g_rp_subvals.append(merged_g_absmag-merged_rp_absmag)
                m_g_subvals.append(merged_g_absmag)
                plt.plot(g_rp_subvals, m_g_subvals, color='magenta', linestyle= None, marker='o', markersize= 6)
        #plt.plot(g_rp_vals, m_g_vals, color='magenta', linestyle= None, marker='o', markersize= 6)
        return g_rp_vals, m_g_vals

    return

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



def plot_bkg_cmd(generic_table= generic_table, absmag='g', colours=['bp','rp']):
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
    polything = plt.hexbin(generic_colour_dif, generic_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
    counts = polything.get_array()
    print(counts.shape)
    counts= np.sqrt(counts)
    polything.set_array(counts)
    polything.autoscale()
    plt.ylim(axes_y)
    plt.gca().invert_yaxis()
    #plt.xlabel(r'$G_{BP} - G_{RP}$')
    plt.xlabel(colours[0]+"-"+colours[1])
    plt.xlim(axes_x)
    #plt.ylabel(r'$M_G$')
    plt.ylabel('M_'+absmag)
    plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)
    return


#def make_cmd(target_table=target_table, generic_table= generic_table, absmag='g', colours=['bp', 'rp']):
    #plot_bkg_cmd(generic_table=generic_table, absmag=absmag, colours=colours)
    #plot_target_table(target_table, absmag=absmag, colours=colours)
    #plt.show()
    #return
#plot_bkg_cmd(colours=['g','rp'])
#merge_stars()
#plt.show()



#merge_stars(star1_table= input_table[0], star2_table=input_table[1],real_stars=False)
##plot_target_table(input_table, colours= ['g','rp'])
#plot_bkg_cmd(colours=['g','rp'])
#plt.show()

#merge_stars(star1_table= input_table[0], star2_table=input_table[1],real_stars=True)
#plot_target_table(input_table, colours= ['g','rp'])
#plot_bkg_cmd(colours=['g','rp'])
#plt.show()
 
star2_dict= find_star2(star1_table, merged_star_table, halfreal=True)
#plot_target_table(input_table, colours= ['g','rp'])
plot_bkg_cmd(colours=['g','rp'])
plt.show()

merge_stars(star1_table=star1_table, star2_table=star2_table, real_stars=True)
plot_bkg_cmd(colours=['g','rp'])
plt.show()

star2_dict= find_star2(star1_table, merged_star_table)
#plot_bkg_cmd()
plot_bkg_cmd(colours=['g','rp'])
plt.show()


star2_dict= find_star2(star2_table, merged_star_table)
plot_bkg_cmd(colours=['g','rp'])
#plot_bkg_cmd()
plt.show()
