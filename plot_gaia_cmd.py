"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-06-11
Take tables output by retrieve_data.py and plot them as the color-magnitude diagram of the gaia data

"""
import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, QTable
import matplotlib.pyplot as plt
import scipy.stats as scistats
import seaborn as sns

num_targs = 'all'
distance = 50
grid_num = 220
mc_number = 10000
percent_off = 34 #1-sigma equivalent
target_label = "PSRJ1435-4715"
#other_target_label= "PSRJ1816+4510"
other_target_label= "PSRJ1023+0038"

if num_targs == 'all':
    print('Distance-limited Sample like Figure 6 from DR2HRD')
    generic_input = 'all_'+str(distance)+'pc_gaia.csv'
    title_suffix = ' in the ' + str(distance)+ 'pc sample following DR2HRD Figure 6'
else:
    num_targs = int(num_targs)
    #generic_input = 'top'+str(num_targs) + '_nearby_gaia.csv'
    #title_suffix= str(num_targs)+ 'stars in the'  +str(distance) + 'pc Gaia CMD'
    title_suffix = 'in the ' +str(num_targs)+ ' star sample following DR2HRD Figure 1'
    generic_input = 'top'+str(num_targs) + '_' +str(distance)+'pc_gaia.csv'
target_input = 'target_gaia.csv'
#other_target_input = "PSRJ1816p4510_gaia.csv"
other_target_input = "PSRJ1023p0038_gaia.csv"


#generic_table = Table.read('top500_nearby_gaia.csv')
#generic_table = Table.read('top5000_nearby_gaia.csv')
generic_table = Table.read(generic_input)
target_table = Table.read(target_input)
other_target_table = Table.read(other_target_input)


def distance_modulus(g_mag, distance, extinction = 0.0):
    return g_mag - 5*np.log10(distance/10.)
    #return g_mag - 5*np.log10(distance/10.)- np.float_(extinction)

def get_mc_distribution(value, error):
    error_distribution = np.random.normal(loc= value, scale = error, size = mc_number)
    return error_distribution
#generic_table.pprint()

def remove_negative(array):
    output_array = array[np.where(array>0)]
    print('Removed ' +str(array.shape[0]-output_array.shape[0]) + ' negatives')
    return output_array

def get_errors(distribution, percent_off = percent_off):
    """
    values for the error bars on the plot
    
    Returns
    [+ value, - value], so to get the points on the plot where they'd be located you do
    value + [+ value] , value - [- value]
    Basically these are the width of the uncertainty range on either side.
    """
    low_bar = np.nanpercentile(distribution, 50-percent_off)
    median = np.nanmedian(distribution)
    high_bar = np.nanpercentile(distribution, 50+percent_off)
    return np.array([[median-low_bar],[high_bar-median]])

try:
    generic_parallax = generic_table ['parallax']
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

print(generic_g_absmag.shape)
#print(target_table['ra'])
#print(target_table['dec'])
target_parallax = target_table['parallax']
target_parallax = target_parallax*1e-3
target_distance = 1./target_parallax
target_parallax_error = target_table['parallax_error']*1e-3
target_parallax_dist = get_mc_distribution(target_parallax, target_parallax_error)
target_parallax_dist = remove_negative(target_parallax_dist)
#remove_negative(target_parallax_dist)
target_distance_dist = 1./target_parallax_dist

#target_distance_err_bounds = 

target_extinction= target_table['a_g_val']

target_g_mag = target_table['phot_g_mean_mag']
target_bp_rp = target_table['bp_rp']





#generic_g_absmag = distance_modulus(generic_g_mag, generic_distance, extinction = generic_extinction)
target_g_absmag = distance_modulus(target_g_mag, target_distance,extinction= target_extinction)
target_g_absmag_dist = distance_modulus(target_g_mag, target_distance_dist, extinction = target_extinction)
target_g_absmag_err = get_errors(target_g_absmag_dist)

plt.hist(target_g_absmag_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
plt.axvline(np.nanmedian(target_g_absmag_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
plt.axvline(np.nanpercentile(target_g_absmag_dist, 84), color = 'cyan')
plt.errorbar( target_g_absmag, 0.5, xerr = target_g_absmag_err, marker = '*', markersize = 8, color = 'b', label = target_label, capsize = 4)
#plt.axvline(x=target_g_absmag, color = 'r', linestyle = ':', label = 'Measured value')
plt.xlabel(r'$M_G$')
plt.title(target_label)
plt.legend()
plt.show()


#print(generic_parallax)
print(target_g_mag)
print(target_extinction)
print(target_table['a_g_percentile_lower'],target_table['a_g_percentile_upper'])
print(target_g_absmag)
print(target_bp_rp)
print(target_g_absmag_err)

########3 other



other_target_parallax = other_target_table['parallax']
other_target_parallax = other_target_parallax*1e-3
other_target_distance = 1./other_target_parallax
other_target_parallax_error = other_target_table['parallax_error']*1e-3
other_target_parallax_dist = get_mc_distribution(other_target_parallax, other_target_parallax_error)
other_target_parallax_dist = remove_negative(other_target_parallax_dist)
#remove_negative(target_parallax_dist)
other_target_distance_dist = 1./other_target_parallax_dist

#target_distance_err_bounds = 

other_target_extinction= other_target_table['a_g_val']

other_target_g_mag =other_target_table['phot_g_mean_mag']
other_target_bp_rp = other_target_table['bp_rp']





#generic_g_absmag = distance_modulus(generic_g_mag, generic_distance, extinction = generic_extinction)
other_target_g_absmag = distance_modulus(other_target_g_mag, other_target_distance,extinction= other_target_extinction)
other_target_g_absmag_dist = distance_modulus(other_target_g_mag, other_target_distance_dist, extinction = other_target_extinction)
other_target_g_absmag_err = get_errors(other_target_g_absmag_dist)

plt.hist(other_target_g_absmag_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
plt.axvline(np.nanmedian(other_target_g_absmag_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
plt.errorbar( other_target_g_absmag, 0.5, xerr = other_target_g_absmag_err, marker = '*', markersize = 8, color = 'b', label = 'Measured value', capsize = 4, linestyle = 'none')
#plt.axvline(x=target_g_absmag, color = 'r', linestyle = ':', label = 'Measured value')
plt.xlabel(r'$M_G$')
plt.title(other_target_label)
plt.legend()
plt.show()


####### 

def make_density_plot(g_abs, bp_rp):
    #Calculate the point density
    xy = np.vstack([np.array(bp_rp),np.array(g_abs)])
    print('starting KDE')
    z = scistats.gaussian_kde(xy)(xy)
    print('finished KDE')
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    g_abs, bp_rp, z = g_abs[idx], bp_rp[idx], z[idx]
    z= np.sqrt(z)
    plt.scatter(bp_rp, g_abs, c=z, s=4, edgecolor = '', cmap= 'hot')
    #fig, ax = plt.subplots()
    #ax.scatter(x, y, c=z, s=50, edgecolor='')
    #plt.show()
    #plt.show()

#plt.scatter(generic_bp_rp, generic_g_absmag, s = 1, alpha = 0.05, color = 'k')
make_density_plot(generic_g_absmag, generic_bp_rp)
#plt.plot(target_bp_rp, target_g_absmag, marker = '*', markersize = 8, color = 'b')
plt.errorbar(target_bp_rp, target_g_absmag, yerr = target_g_absmag_err, marker = '*', markersize = 8, color = 'b', capsize = 4, label = target_label, linestyle = 'none')
plt.errorbar(other_target_bp_rp, other_target_g_absmag, yerr = other_target_g_absmag_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')
#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
#counts = polything.get_array()
#print(counts.shape)
#counts= np.sqrt(counts)
#polything.set_array(counts)
#polything.autoscale()
#plt.hist2d(generic_bp_rp, generic_g_absmag, bins = 100, cmap = 'Reds')
#plt.title('PSR J1431-4715' + title_suffix)
plt.title('PSR J1431-4715' + title_suffix)

#plt.ylim([-4, 16])

plt.gca().invert_yaxis()
plt.xlabel(r'$G_{BP} - G_{RP}$')
#plt.xlim([-1,5])
plt.ylabel(r'$M_G$')
plt.legend()
plt.show()

#plt.hist2d(generic_bp_rp, generic_g_absmag, bins = 100, cmap = 'Reds')
#plt.gca().invert_yaxis()
#plt.xlabel(r'$G_{BP} - G_{RP}$')
#plt.ylabel(r'$G$')
#plt.show()


#plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
#plt.gca().invert_yaxis()
#plt.xlabel(r'$G_{BP} - G_{RP}$')
#plt.ylabel(r'$G$')
#plt.show()

#sns.kdeplot(np.array(generic_bp_rp), np.array(generic_g_absmag), cmap = 'Reds', shade = True)
#plt.show()
