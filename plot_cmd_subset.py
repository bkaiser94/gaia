"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-09-19

Basically the same idea as plot_gaia_cmd.py, but this one should be plotting the primed coordinates from the transform to the new basis



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


input_filename= 'test_subset200.txt'
output_name= 'subset_weirdos.txt'
other_output_name='20190107_chris.csv'

generic_table = Table.read(input_filename, format= 'ascii.tab')

axes_x= [-0.7, 5]
#axes_x= [-0.7, 2]
#axes_x= [-2, 5]


axes_y = [-2,16]
#axes_y = [-4,18]
#axes_y = [9,12]
#axes_y= [9,16]

grid_num=1000
percent_off = 34 #1-sigma equivalent
mc_number = 10000


basis_coeff1= 0.968664
basis_coeff2= 0.248375

bp_rp_bounds = [1.54,2.40]
gabsmag_bounds=[13.78,14.6]


def distance_modulus(g_mag, distance, extinction = 0.0):
    return g_mag - 5*np.log10(distance/10.)
    #return g_mag - 5*np.log10(distance/10.)- np.float_(extinction)
    
def get_mag(flux, filter_string):
    mag0 = zeropoint_dict[filter_string][0]
    return -2.5*np.log10(flux) +mag0

def get_mc_distribution(value, error):
    error_distribution = np.random.normal(loc= value, scale = error, size = mc_number)
    return error_distribution
#generic_table.pprint()

def remove_negative(array, verbose= True):
    output_array = array[np.where(array>0)]
    if (verbose and array.shape[0]-output_array.shape[0] >0):
        print('Removed ' +str(array.shape[0]-output_array.shape[0]) + ' negatives')
    return output_array

def match_sizes(change_array, match_array):
    """
    Intended to keep compatibility with an array that has had negatives removed
    """
    #if ((type(change_array) == type(np.array([0]))) and (type(match_array) == type(np.array([0])))):
        #print("type match")
    try:
        min_inds = np.nanmin([change_array.shape[0], match_array.shape[0]])
        
        return change_array[:min_inds], match_array[:min_inds]
    except AttributeError:
        #the inputs aren't actually arrays
        return change_array, match_array
    #else:
        #return change_array, match_array
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
        #print(error)
        return  np.array([[median.value-low_bar],[high_bar-median.value]])
    #return np.vstack([np.array(median-low_bar),np.array(high_bar-median)])


def get_filter_vals(table, filter_string):
    flux_string = 'phot_'+filter_string+ '_mean_flux'
    #print(flux_string)
    #print("table", type(table))
    #print(table)
    phot_mean_flux = table[flux_string]
    error_string = flux_string + '_error'
    phot_mean_flux_error = table[error_string]
    flux_distribution = get_mc_distribution(phot_mean_flux, phot_mean_flux_error)
    return phot_mean_flux, flux_distribution



def get_bp_rp(table, plot_all = False, verbose =True):
    bp_mean_flux, bp_dist = get_filter_vals(table, 'bp')
    rp_mean_flux, rp_dist = get_filter_vals(table, 'rp')
    bp_mag = get_mag(bp_mean_flux, 'bp')
    rp_mag = get_mag(rp_mean_flux, 'rp')
    if verbose:
        print("bp_calc-bp_measured", bp_mag - table['phot_bp_mean_mag'])
        print("rp_calc - rp_measured", rp_mag - table['phot_rp_mean_mag'])
    bp_mag_dist = get_mag(bp_dist, 'bp')
    rp_mag_dist = get_mag(rp_dist, 'rp')
    bp_rp = bp_mag- rp_mag
    bp_rp_dist= bp_mag_dist- rp_mag_dist
    bp_rp_error = get_errors(bp_rp_dist)
    if plot_all:
        plt.hist(bp_rp_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
        plt.axvline(np.nanmedian(bp_rp_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
        plt.axvline(np.nanpercentile(bp_rp_dist, 84), color = 'cyan')
        plt.errorbar(bp_rp, 0.5, xerr = bp_rp_error, marker = '*', markersize = 8, color = 'b', label = "BP-RP", capsize = 4)
        #plt.axvline(x=target_g_absmag, color = 'r', linestyle = ':', label = 'Measured value')
        plt.xlabel(r'$G_{BP}-G_{RP}$')
        #plt.title(target_label)
        plt.legend()
        plt.show()
    return bp_rp, bp_rp_error

#def slice_and_dice(array1, array2, bounds1, bounds2):
    #lower_cut= np.where(array1 > bounds1[0])
    #low_array1 = array1[lower_cut]
    #low_array2= array2[lower_cut]
    #upper_cut = np.where(low_array1< bounds1[1])
    #upper_array1= low_array1[upper_cut]
    #upper_array2= low_array2[upper_cut]
    #left_cut= np.where(upper_array2 > bounds2[0])
    #left_array1= upper_array1[left_cut]
    #left_array2= upper_array2[left_cut]
    #right_cut= np.where(left_array2 < bounds2[1])
    #right_array1= left_array1[right_cut]
    #right_array2= left_array2[right_cut]
    #return right_array1, right_array2
    
def slice_and_dice(table, bounds1, bounds2):
    lower_cut= np.where(table['xprime'] > bounds1[0])
    low_table= table[lower_cut]
    upper_cut = np.where(low_table['xprime'] < bounds1[1])
    upper_table= low_table[upper_cut]
    left_cut= np.where(upper_table['yprime'] > bounds2[0])
    left_table= upper_table[left_cut]
    right_cut= np.where(left_table['yprime'] < bounds2[1])
    right_table= left_table[right_cut]
    return right_table


def chop_chop(table, bounds1, bounds2):
    lower_cut= np.where(table['bp_rp'] > bounds1[0])
    low_table= table[lower_cut]
    upper_cut = np.where(low_table['bp_rp'] < bounds1[1])
    upper_table= low_table[upper_cut]
    left_cut= np.where(upper_table['mg'] > bounds2[0])
    left_table= upper_table[left_cut]
    right_cut= np.where(left_table['mg'] < bounds2[1])
    right_table= left_table[right_cut]
    return right_table
   
try:
    generic_parallax = generic_table ['parallax']+parallax_correction
    generic_parallax = generic_parallax *1e-3 #parallax in arcseconds now
    generic_distance = 1./generic_parallax #parsec distance

    generic_extinction = generic_table['a_g_val']

    generic_g_mag = generic_table['phot_g_mean_mag']
    generic_bp_rp = generic_table['bp_rp']
    generic_g_absmag = distance_modulus(generic_g_mag, generic_distance, extinction = generic_extinction)

except KeyError as error:
    print(error)
    print("assuming it's the simplified file.")
    generic_g_absmag= generic_table['mg']
    generic_bp_rp= generic_table['bp_rp']
    generic_xprime= generic_table['xprime']
    generic_yprime= generic_table['yprime']
    
#print(generic_table[0])
#in_bounds = np.where(generic_xprime> 0.877)
acc_table = slice_and_dice(generic_table,[0.877,2.278], [4.627,9])
print(acc_table)


acc_table.write(output_name, format = 'ascii.tab', overwrite= True)



def plot_context(acc_table= acc_table):
    polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
    counts = polything.get_array()
    #print(counts.shape)
    counts= np.sqrt(counts)
    polything.set_array(counts)
    polything.autoscale()

    plt.plot([0.3,0.3-basis_coeff1],[4.732, 4.732+basis_coeff2])
    plt.plot([0.3,0.3+basis_coeff2],[4.732, 4.732+basis_coeff1])
    #plt.scatter(generic_bp_rp, generic_g_absmag)

    plt.scatter(acc_table['bp_rp'], acc_table['mg'], color='b')

    #plt.title('PSR J1431-4715' + title_suffix)

    #plt.ylim([-4, 16])
    plt.ylim(axes_y)
    #plt.ylim([-4,16])

    plt.gca().invert_yaxis()
    plt.xlabel(r'$G_{BP} - G_{RP}$')
    #plt.xlim([-1,5])
    plt.xlim(axes_x)
    #plt.xlim([-4,16])
    plt.ylabel(r'$M_G$')
    #plt.legend()
    plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)


################ The search region ###################
middle_area= generic_table[np.where(generic_table['bp_rp']>2.25)]
middler_area= middle_area[np.where(middle_area['bp_rp']<3.7)]
only_one= middler_area[np.where(middler_area['mg']> 14.5)]
print("hopefully low point")
print(only_one)

#chris_table= chop_chop(generic_table, [2.25, 3.7],[14.5, 20])
chris_table= chop_chop(generic_table, bp_rp_bounds, gabsmag_bounds)
chris_sort = np.argsort(chris_table['bp_rp'])
chris_table= chris_table[chris_sort]
chris_table.write(other_output_name, format='csv')
for row in chris_table:
    print(row)
    plot_context(acc_table= chris_table)
    plt.scatter(row['bp_rp'], row['mg'], color= 'r')
    plt.show()
############################################

param_distance = acc_table['xprime']**2+acc_table['yprime']**2
sort_order = np.argsort(param_distance)
sorted_table= acc_table[sort_order]
#for row in acc_table:
    #print(row)
    #plot_context()
    #plt.scatter(row['bp_rp'], row['mg'], color= 'r')
    #plt.show()

for row in sorted_table:
    print(row)
    plot_context()
    plt.scatter(row['bp_rp'], row['mg'], color= 'r')
    plt.show()


polything = plt.hexbin(generic_xprime, generic_yprime, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
counts = polything.get_array()
#print(counts.shape)
counts= np.sqrt(counts)
polything.set_array(counts)
polything.autoscale()

#plt.scatter(generic_xprime, generic_yprime)
plt.scatter(acc_table['xprime'], acc_table['yprime'], color='b')

#plt.title('PSR J1431-4715' + title_suffix)

#plt.ylim([-4, 16])
#plt.ylim(axes_y)

#plt.gca().invert_yaxis()
plt.xlabel(r"$x'$")
plt.axvline(x=0.877)
plt.axvline(x=2.278)
plt.axhline(y=4.627)
plt.axhline(y=9)
#plt.xlim([-1,5])
#plt.xlim(axes_x)
plt.ylabel(r"$y'$")
#plt.legend()
plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)
plt.show()
