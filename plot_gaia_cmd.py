"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-06-11
Take tables output by retrieve_data.py and plot them as the color-magnitude diagram of the gaia data

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

import passband_model_convolution as pmc
import wdatmos

plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
#plt.rc('font', size = 11)
plt.rc('lines', markersize = 5)
precision = 2
parallax_correction = 0.029 #from Lindgren et al 2018
#parallax_correction = 0 #so nothing done

teff = 7250
logg = 6.0

q= 10.5823898006 #from RV fitting (it's pretty loose)
#test_radii = np.array([0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 100]) * u.Rsun
test_radii = np.array([0.01, 0.1, 1, 10]) * u.Rsun

num_targs = 'all'
distance = 100
#distance = 25
#grid_num = 220
grid_num = 225
mc_number = 10000
percent_off = 34 #1-sigma equivalent


target_label = "PSRJ1435-4715"
#other_target_label= "PSRJ1816+4510"
#other_target_label= "PSRJ1023+0038"
#other_target_label = "Crab Pulsar"
other_target_label = "PSRJ1435-6100"


target_input = 'target_gaia.csv'
#other_target_input = "PSRJ1816p4510_gaia.csv"
#other_target_input = "PSRJ1023p0038_gaia.csv"
#other_target_input = "Crab_gaia.csv"
other_target_input= "PSRJ1435m6100.csv"
#other_target_input= 'PSRJ1903p0327.csv'


zeropoint_dict={"g": [25.6883657251, 0.0017850023],
                "bp": [ 25.3513881707 , 0.0013918258],
                "rp": [24.7619199882, 0.0019145719]} #from Evans et al 2018, the DR2 values [ZP, sigma]




#zeropoint_dict={"g": [ 25.7933969562,  0.0017848281],
                #"bp": [25.3805596387,  0.0013917453], 
                #"rp": [25.1161276701, 0.001914645] } #AB from Evans et al 2018, the DR2 values [ZP, sigma]
                
#extinction coefficients from DR2HRD 2018 Table 1
#the 0th entry is blank to match the indexing of the table
c_coeffs= {'g': [0, 0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0094],
           'bp': [0, 1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043],
           'rp': [0, 0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005,0.0006]}
                

dtype_list = ['S32', 'float', 'float', 'float', 'float','float', 'float', 'float', 'float', 'float','float', 'float', 'float', 'float', 'float']
pulsar_list_file = 'Jennings_table2.txt'
pulsar_list_all = np.genfromtxt(pulsar_list_file, delimiter = '\t', names = True, dtype = dtype_list)



if num_targs == 'all':
    print('Distance-limited Sample like Figure 6 from DR2HRD')
    #generic_input = 'all_'+str(distance)+'pc_gaia.csv'
    generic_input = 'all_'+str(distance)+'pc_gaia_corr.csv'
    #title_suffix =' and ' + other_target_label+ '(green) in ' + str(distance)+ 'pc (DR2HRD Figure 6)'
    #title_suffix = ' and ' + other_target_label+ '(green) in ' + str(distance)+ 'pc (DR2HRD Figure 6) (0.029 mas correction)''
    title_suffix = ' in the ' + str(distance)+ 'pc Gaia DR2 CMD'

else:
    num_targs = int(num_targs)
    generic_input = 'top'+str(num_targs) + '_nearby_gaia.csv'
    #title_suffix= str(num_targs)+ 'stars in the'  +str(distance) + 'pc Gaia CMD'
    title_suffix = 'in the ' +str(num_targs)+ ' star sample following DR2HRD Figure 1'
    #generic_input = 'top'+str(num_targs) + '_' +str(distance)+'pc_gaia.csv'



#generic_table = Table.read('top500_nearby_gaia.csv')
#generic_table = Table.read('top5000_nearby_gaia.csv')
generic_table = Table.read(generic_input)
target_table = Table.read(target_input)
other_target_table = Table.read(other_target_input)


def distance_modulus(g_mag, distance, extinction = 0.0):
    return g_mag - 5*np.log10(distance/10.)
    #return g_mag - 5*np.log10(distance/10.)- np.float_(extinction)
    
def get_mag(flux, filter_string):
    mag0 = zeropoint_dict[filter_string][0]
    return -2.5*np.log10(flux) +mag0

def get_mc_distribution(value, error):
    error_distribution = np.random.normal(loc= value, scale = error, size = mc_number)
    #try:
        #print(value.shape)
        #print(error.shape)
        #error_distribution = np.random.normal(loc= value, scale = error, size = (value.shape[0], mc_number))
    #except TypeError:
        #error_distribution = np.random.normal(loc= value, scale = error, size = mc_number)
    return error_distribution
#generic_table.pprint()

def remove_negative(array, verbose= True):
    output_array = array[np.where(array>0)]
    if (verbose and array.shape[0]-output_array.shape[0] >0):
        print('Removed ' +str(array.shape[0]-output_array.shape[0]) + ' negatives')
    return output_array

def get_errors(distribution, percent_off = percent_off):
    """
    values for the error bars on the plot
    
    Returns
    [+ value, - value], so to get the points on the plot where they'd be located you do
    value + [+ value] , value - [- value]
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


def get_a_x(bp_rp,  a_0, passband_string = 'g'):
    """
    Equation 1 from DR2HRD 2018
    """
    c_vals = c_coeffs[passband_string]
    k_x=1
    k_x= c_vals[1]+c_vals[2]*bp_rp+c_vals[3]*bp_rp**2+c_vals[4]*bp_rp**3+c_vals[5]*a_0+c_vals[6]*a_0**2+c[7]*bp_rp*a_0
    a_x = a_0*k_x
    return a_x

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



def get_g_abs_mag(table, plot_all = False):
    g_mean_flux, g_dist = get_filter_vals(table, 'g')
    g_mag = get_mag(g_mean_flux, 'g')
    g_mag_dist = get_mag(g_dist, 'g')
    print("g_calc-g_measured", g_mag - table['phot_g_mean_mag'])
    parallax = table['parallax']+parallax_correction
    parallax = parallax*1e-3
    distance = 1./parallax
    parallax_error = table['parallax_error']*1e-3
    parallax_dist = get_mc_distribution(parallax, parallax_error)
    parallax_dist = remove_negative(parallax_dist)
    if parallax < 0:
        parallax_median = np.nanmedian(parallax_dist)
        print("PARALLAX < 0!", parallax, "setting to median of positive distribution:", parallax_median)
        parallax = parallax_median
    distance = 1./parallax
    distance_dist = 1./parallax_dist
    index_length = distance_dist.shape[0]
    #print("g_mag_dist.shape", g_mag_dist.shape, "distance_dist.shape", distance_dist.shape)
    g_mag_dist = g_mag_dist[:index_length]
    #print("g_mag_dist.shape", g_mag_dist.shape, "distance_dist.shape", distance_dist.shape)
    g_abs_mag = distance_modulus(g_mag, distance)
    g_abs_mag_dist = distance_modulus(g_mag_dist, distance_dist)
    g_abs_mag_error= get_errors(g_abs_mag_dist)

    #extinction= table['a_g_val']

    if plot_all:
        plt.hist(g_abs_mag_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
        plt.axvline(np.nanmedian(g_abs_mag_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
        plt.axvline(np.nanpercentile(g_abs_mag_dist, 84), color = 'cyan')
        plt.errorbar(g_abs_mag, 0.5, xerr = g_abs_mag_error, marker = '*', markersize = 8, color = 'b', label = r"$M_G$", capsize = 4)
        #plt.axvline(x=target_g_absmag, color = 'r', linestyle = ':', label = 'Measured value')
        plt.xlabel(r'$M_{G}$')
        #plt.title(target_label)
        plt.legend()
        plt.show()
    return g_abs_mag, g_abs_mag_error, g_abs_mag_dist

def correct_bp_flux(table, teff= teff, logg=logg):
    obs_bp_flux, obs_bp_flux_dist = get_filter_vals(table, 'bp')
    obs_rp_flux, obs_rp_flux_dist = get_filter_vals(table, 'rp')
    
    return


def calc_ns_mass(comp_mass, q=q):
    return q*comp_mass

def calc_mean_density(mass, radius):
    return (mass/(4./3 *np.pi*radius**3)).cgs
#def get_mass(radius, logg=logg):
    #logg = logg*u.
    
def calc_min_period(mass, radius):
    """
    Returns the minimum period in hours that a given density star could have. 
    immediately converted to days though, so it actually returns the value in days
    """
    mean_density = calc_mean_density(mass, radius).cgs.value
    return ((107./mean_density)**2*u.hour).to(u.day)

def get_pass_abs_mag(table, plot_all = False, passband_string= 'g', verbose = True):
    mean_flux, flux_dist = get_filter_vals(table, passband_string)
    mag = get_mag(mean_flux, passband_string)
    mag_dist = get_mag(mean_flux, passband_string )
    
    parallax = table['parallax']+parallax_correction
    parallax = parallax*1e-3
    distance = 1./parallax
    parallax_error = table['parallax_error']*1e-3
    parallax_dist = get_mc_distribution(parallax, parallax_error)
    parallax_dist = remove_negative(parallax_dist, verbose= verbose)
    if verbose:
        print(passband_string+ "_calc" + "-" + passband_string+ "_measured", mag - table['phot_' +passband_string+'_mean_mag'])
    if parallax < 0:
        parallax_median = np.nanmedian(parallax_dist)
        if verbose:
            print("PARALLAX < 0!", parallax, "setting to median of positive distribution:", parallax_median)
        parallax = parallax_median
    distance = 1./parallax
    distance_dist = 1./parallax_dist
    index_length = distance_dist.shape[0]
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
    return abs_mag, abs_mag_error, abs_mag_dist

def get_rad_mass(table, teff=teff, logg=logg, passband_string= 'g', plot_all = False, verbose= True):
    """
    Get the radius and mass for a given model when provided with the given absolute magnitude
    for a given band. Return the radius and mass for that absolute magnitude (or distribution) for the
    provided model.
    """
    abs_mag,abs_mag_err, abs_mag_dist = get_pass_abs_mag(table, plot_all = plot_all, passband_string=passband_string, verbose = verbose)
    radius = pmc.get_radius(abs_mag, teff= teff, logg= logg, passband_string= passband_string)
    mass = (pmc.get_mass(radius, logg)).to(u.Msun)
    radius_dist = pmc.get_radius(abs_mag_dist, teff = teff, logg=logg, passband_string = passband_string)
    mass_dist = (pmc.get_mass(radius_dist, logg)).to(u.Msun)
    radius_err = get_errors(radius_dist)
    mass_err= get_errors(mass_dist)
    if plot_all:
        plt.hist(radius_dist.value, bins=75, normed=1, label = 'MC Distribution', color = 'g')
        plt.axvline(np.nanmedian(radius_dist.value), color = 'r', linestyle = '--', label = 'Median of MC Dist')
        plt.axvline(np.nanpercentile(radius_dist.value, 84), color = 'cyan')
        plt.errorbar(radius.value, 0.5, xerr = radius_err, marker = '*', markersize = 8, color = 'b', label = r"$R_{*}$", capsize = 4)
        #plt.axvline(x=target_g_absmag, color = 'r', linestyle = ':', label = 'Measured value')
        plt.xlabel(r'$R_{*}(R_{\odot})$')
        plt.title("R calculated from absolute " + passband_string + " magnitude")
        plt.legend()
        plt.show()
        plt.hist(mass_dist.value, bins=1000, label = 'MC Distribution', color = 'g')
        plt.axvline(np.nanmedian(mass_dist.value), color = 'r', linestyle = '--', label = 'Median of MC Dist')
        plt.axvline(np.nanpercentile(mass_dist.value, 84), color = 'cyan')
        plt.errorbar(mass.value, 100, xerr = mass_err, marker = '*', markersize = 8, color = 'b', label = r"$M_{*}$", capsize = 4)
        #plt.axvline(x=target_g_absmag, color = 'r', linestyle = ':', label = 'Measured value')
        plt.xlabel(r'$M_{*}(M_{\odot})$')
        plt.title(r"$M_{*}$ calculated from absolute " + passband_string + " magnitude")
        plt.legend()
        plt.show()
    mean_density = calc_mean_density(mass, radius)
    mean_density_dist = calc_mean_density(mass_dist, radius_dist)
    mean_density_err = get_errors(mean_density_dist)
    ns_mass = calc_ns_mass(mass)
    ns_mass_dist= calc_ns_mass(mass_dist)
    ns_mass_err= get_errors(ns_mass_dist)
    min_period = calc_min_period(mass, radius)
    min_period_dist= calc_min_period(mass_dist, radius_dist)
    min_period_err = get_errors(min_period_dist)
    print(passband_string+ " Abs Mag:", abs_mag[0], "-/+", abs_mag_err[0,0], abs_mag_err[1,0])
    print(passband_string+ " Radius:", radius[0], "-/+", radius_err[0,0], radius_err[1,0])
    print(passband_string+ " Comp Mass:", mass[0], "-/+", mass_err[0,0],mass_err[1,0])
    print(passband_string+ " PSR Mass:",ns_mass[0], "-/+", ns_mass_err[0,0],ns_mass_err[1,0])
    print(passband_string+ " Mean Density:", mean_density[0], "-/+", mean_density_err[0,0], mean_density_err[1,0])
    print(passband_string+ " Min Period:", min_period[0], "-/+", min_period_err[0,0], min_period_err[1,0])

    return radius, mass, radius_dist, mass_dist



    
    
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



#target_g_absmag, target_g_absmag_err, target_g_absmag_dist = get_g_abs_mag(target_table, plot_all = True)
target_g_absmag, target_g_absmag_err, target_g_absmag_dist = get_pass_abs_mag(target_table, plot_all = True, passband_string = 'g')
target_bp_rp, target_bp_rp_err= get_bp_rp(target_table, plot_all = True)

#target_g_radius = pmc.get_radius(target_g_absmag, teff= teff, logg= logg, passband_string= 'G')
#target_g_radius_dist = pmc.get_radius(target_g_absmag_dist, teff= teff, logg= logg, passband_string= 'G')


#target_g_radius_err = get_errors(target_g_radius_dist)

#target_g_mass = (pmc.get_mass(target_g_radius, logg)).to(u.Msun)
#target_g_mass_dist= (pmc.get_mass(target_g_radius_dist, logg)).to(u.Msun)
#target_g_mass_err = get_errors(target_g_mass_dist)

#target_g_radius, target_g_mass = get_rad_mass(target_g_absmag, logg=logg, teff= teff, passband_string= 'g')
target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all = True)
#target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_g_absmag_dist, logg=logg, teff= teff, passband_string= 'g')
target_g_radius_err = get_errors(target_g_radius_dist)
target_g_mass_err = get_errors(target_g_mass_dist)


sim_target_gabsmag, sim_target_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff= teff, radius =target_g_radius)
sim_target_gabsmag_dist, sim_target_bp_rp= pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = target_g_radius_dist)
sim_target_gabsmag_err = get_errors(sim_target_gabsmag_dist)


print("Target Radius:", target_g_radius, "-/+", target_g_radius_err)
print("Target Mass:", target_g_mass, "-/+", target_g_mass_err)



model_g_absmag, model_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = test_radii)


########3 other

###### BP

target_bp_radius, target_bp_mass,target_bp_radius_dist, target_bp_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'bp', plot_all = True)
#target_bp_radius_dist, target_bp_mass_dist = get_rad_mass(, logg=logg, teff= teff, passband_string= 'bp')
target_bp_radius_err = get_errors(target_bp_radius_dist)
target_bp_mass_err = get_errors(target_bp_mass_dist)

print("Target Radius bp:", target_bp_radius, "-/+", target_bp_radius_err)
print("Target Mass bp:", target_bp_mass, "-/+", target_bp_mass_err)


sim_target_bp_absmag, sim_target_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff= teff, radius =target_bp_radius)
sim_target_bp_absmag_dist, sim_target_bp_rp= pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = target_bp_radius_dist)
sim_target_bp_absmag_err = get_errors(sim_target_bp_absmag_dist)

#####3 RP



target_rp_radius, target_rp_mass,target_rp_radius_dist, target_rp_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'rp', plot_all = True)
#target_rp_radius_dist, target_rp_mass_dist = get_rad_mass(target_rp_absmag_dist, logg=logg, teff= teff, passband_string= 'rp')
target_rp_radius_err = get_errors(target_rp_radius_dist)
target_rp_mass_err = get_errors(target_rp_mass_dist)

print("Target Radius rp:", target_rp_radius, "-/+", target_rp_radius_err)
print("Target Mass rp:", target_rp_mass, "-/+", target_rp_mass_err)


sim_target_rp_absmag, sim_target_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff= teff, radius =target_rp_radius)
sim_target_rp_absmag_dist, sim_target_bp_rp= pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = target_rp_radius_dist)
sim_target_rp_absmag_err = get_errors(sim_target_rp_absmag_dist)



################
#other_target_g_absmag, other_target_g_absmag_err, other_target_g_absmag_dist= get_g_abs_mag(other_target_table, plot_all = True)
#other_target_bp_rp, other_target_bp_rp_err = get_bp_rp(other_target_table, plot_all = True)
###################


#plt.hist(other_target_g_absmag_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
#plt.axvline(np.nanmedian(other_target_g_absmag_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
#plt.errorbar( other_target_g_absmag, 0.5, xerr = other_target_g_absmag_err, marker = '*', markersize = 8, color = 'b', label = 'Measured value', capsize = 4, linestyle = 'none')
##plt.axvline(x=target_g_absmag, color = 'r', linestyle = ':', label = 'Measured value')
#plt.xlabel(r'$M_G$')
#plt.title(other_target_label)
#plt.legend()
#plt.show()


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


plt.errorbar(target_bp_rp, target_g_absmag, yerr = target_g_absmag_err, xerr = target_bp_rp_err, marker = '*', markersize = 8, color = 'b', capsize = 4, label = target_label, linestyle = 'none')

#plt.errorbar(other_target_bp_rp, other_target_g_absmag, yerr = other_target_g_absmag_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')
#This one \/ \/ \/
#plt.errorbar(other_target_bp_rp, other_target_g_absmag, yerr = other_target_g_absmag_err, xerr= other_target_bp_rp_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')

#plt.errorbar(other_target_bp_rp, other_target_g_absmag, xerr= other_target_bp_rp_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')

#plt.plot(model_bp_rp*np.ones(model_g_absmag.shape[0]), model_g_absmag, marker = '*', markersize= 8, color = 'green', label = "Model spectra at various radii " + str(np.round(test_radii.value, precision))+ " logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')

plt.plot(sim_target_bp_rp, sim_target_gabsmag, marker = '*', markersize= 8, color = 'green', label = "R_G= "+str(np.round(target_g_radius, precision)[0]) + " M= "+ str(np.round(target_g_mass, precision)[0])+" logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')

plt.plot(sim_target_bp_rp, sim_target_bp_absmag, marker = '*', markersize= 8, color = 'cyan', label = "R_BP= "+str(np.round(target_bp_radius, precision)[0]) + " M= "+ str(np.round(target_bp_mass, precision)[0])+" logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')

plt.plot(sim_target_bp_rp, sim_target_rp_absmag, marker = '*', markersize= 8, color = 'magenta', label = "R_RP= "+str(np.round(target_rp_radius, precision)[0]) + " M= "+ str(np.round(target_rp_mass, precision)[0])+" logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')

#plt.plot(model_bp_rp*np.ones(model_g_absmag.shape[0]), model_g_absmag, marker = '*', markersize= 8, color = 'cyan', label = "Model spectra at various radii  logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')
#plt.legend()


#plt.plot( mass_bp_rp*np.ones(mass_gabsmag.shape[0]), mass_gabsmag,marker = '*', markersize= 8, color = 'magenta', label = "using the mass ratio and surface gravity" + str(logg) + " Teff: " +str(teff), linestyle = 'none')
plt.legend()

polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
counts = polything.get_array()
print(counts.shape)
counts= np.sqrt(counts)
polything.set_array(counts)
polything.autoscale()

plt.title('PSR J1431-4715' + title_suffix)

plt.ylim([-4, 16])

plt.gca().invert_yaxis()
plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.xlim([-1,5])
plt.ylabel(r'$M_G$')
#plt.legend()
plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)
plt.show()


#model_photo_excess= 

###################################3

wd=wdatmos.wdmodel(filename='ELM.hdf5')
teff_array=wd.Teffs
logg_array = wd.loggs

#for teff,logg in zip(teff_array, logg_array):
for logg,teff in zip( logg_array,teff_array):

    print("Teff:", teff, "logg:", logg)
    model = wd(Teff = teff , logg = logg)
    target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all= False, verbose= False)
    #target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_g_absmag_dist, logg=logg, teff= teff, passband_string= 'g')
    target_g_radius_err = get_errors(target_g_radius_dist)
    target_g_mass_err = get_errors(target_g_mass_dist)
    sim_target_gabsmag, sim_target_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff= teff, radius =target_g_radius)
    sim_target_gabsmag_dist, sim_target_bp_rp= pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = target_g_radius_dist)
    sim_target_gabsmag_err = get_errors(sim_target_gabsmag_dist)
    print( "BP-RP:", sim_target_bp_rp)
    #print("G absmag:", sim_target_gabsmag, "-/+", sim_target_gabsmag_err[0,0], sim_target_gabsmag_err[1,0])
    #print("comp_mass:", target_g_mass, "-/+", get_errors(target_g_mass_dist))
    #print("NS_mass: ", calc_ns_mass(target_g_mass), "-/+", get_errors(calc_ns_mass(target_g_mass_dist)))
    print("--------------------------")
    plt.errorbar(sim_target_bp_rp, sim_target_gabsmag, yerr= sim_target_gabsmag_err, marker = '*', markersize= 8, color = 'green')

        










polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
counts = polything.get_array()
print(counts.shape)
counts= np.sqrt(counts)
polything.set_array(counts)
polything.autoscale()

plt.title('PSR J1431-4715' + title_suffix)

plt.ylim([-4, 16])

plt.gca().invert_yaxis()
plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.xlim([-1,5])
plt.ylabel(r'$M_G$')
#plt.legend()
plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)
plt.show()






