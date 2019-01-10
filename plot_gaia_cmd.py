"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-06-11

Take tables output by retrieve_data.py and plot them as the color-magnitude diagram of the gaia data

I have begun experimenting with 'filled' tables, meaning they aren't masked, but instead have values that can actually be detected. One can also set them to be certain values. It would make sense to make my values be on the edges of the CMD

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
import gaia_extinction
import wdatmos
import plotting_dicts as pod

plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
#plt.rc('font', size = 11)
plt.rc('lines', markersize = 5)



precision = 2

#output_filename= 'model_values.csv'
#output_filename= 'model_values_3,4m.csv'
#output_filename= 'model_values_1,25m.csv'
#output_filename= 'model_values_1,4m.csv'
output_filename= 'model_values_1,4m_photo_koester.csv'

#m_psr = 1.9* u.Msun
#m_psr = 3.4*u.Msun
#m_psr = 1.25*u.Msun
m_psr= 1.4*u.Msun

#ext_method = 'B-V'
#ext_method = 'statistical'
ext_method=  'MG-MRP'

list_color= 'b'
single_list=True#turns off the original for-loop method for plotting from a single list
error_bar=True #turns off error bars on the multiple list plot, meaning it has no effect on anything if single_list==True

parallax_correction = 0.029 #from Lindgren et al 2018
#parallax_correction = 0 #so nothing done


#good_teff_range = [6000, 10000]
good_teff_range = [6000, 14750]


#axes_x= [-0.7, 5]
axes_x=[-2,7]
#axes_x=[-0.4,0.9]
#axes_x=[-4, 8]


#axes_y=[-1, 14.5]
#axes_y = [-2,16]
axes_y = [-5,18]
#axes_y = [-8,20]
#axes_y = [9,12]
#axes_y= [9,16]
g_abs_mag_fill= axes_y[1]
bp_rp_fill= axes_x[0]


#axes_y = [4,8]

#teff = 7250
#logg = 6.0


teff = 7500
logg = 5.0


q= 10.5823898006 #from RV fitting (it's pretty loose)
q_err=0.0531987564067
#test_radii = np.array([0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 100]) * u.Rsun
test_radii = np.array([0.01, 0.1, 1, 10]) * u.Rsun

num_targs = 'all'
#num_targs = '47Tuc'
distance = 500
#distance=500
#distance = 25
#grid_num = 220
grid_num = 225
#grid_num = 500
#grid_num =1000

mc_number = 10000
percent_off = 34 #1-sigma equivalent
#percent_off =  #1-sigma equivalent


target_label = "PSRJ1435-4715"
#target_label = "WD1145"

#other_target_label= "PSRJ1816+4510"
#other_target_label= "PSRJ1023+0038"
#other_target_label = "Crab Pulsar"
#other_target_label = "PSRJ1435-6100"
other_target_label= 'BLAP'


#target_input = 'target_gaia.csv'
target_input= 'weird_ogle.csv'
#target_input= 'WD1145p017.csv'
#other_target_input = "PSRJ1816p4510_gaia.csv"
#other_target_input = "PSRJ1023p0038_gaia.csv"
#other_target_input = "Crab_gaia.csv"
#other_target_input= "PSRJ1435m6100.csv"
#other_target_input= 'BLAPs_gaia.csv'
#other_target_input=  'pulsar_companions_gaia.csv'
#other_target_input='hot_wind_wds_gaia.csv'
#other_target_input='hv_wds_gaia.csv'
other_target_input='pre_elms_gaia.csv'
other_target_input='20190107_chris_merge_gaia.csv'
#other_target_input= 'PSRJ1903p0327.csv'


zeropoint_dict={"g": [25.6883657251, 0.0017850023],
                "bp": [ 25.3513881707 , 0.0013918258],
                "rp": [24.7619199882, 0.0019145719]} #from Evans et al 2018, the DR2 values [ZP, sigma]




#zeropoint_dict={"g": [ 25.7933969562,  0.0017848281],
                #"bp": [25.3805596387,  0.0013917453], 
                #"rp": [25.1161276701, 0.001914645] } #AB from Evans et al 2018, the DR2 values [ZP, sigma]
 
start_header= "teff, logg, corr_g_absmag, corr_g_absmag_err_lo, corr_g_absmag_err_hi, R_c, R_c_err_lo, R_c_err_hi, m_c, m_c_err_lo, m_c_err_hi, m_psr, m_psr_err_lo, m_psr_err_hi, mean_density, mean_density_err_lo, mean_density_err_hi, P_min, P_min_err_lo, P_min_err_hi"

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
elif num_targs== '47Tuc':
    generic_input= "47Tuc_10arcmin.csv"
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

other_target_table= other_target_table.filled()


print("other G")
print(other_target_table['phot_g_mean_mag'])
print('other BP')
print(other_target_table['phot_bp_mean_mag'])
print('other RP')
print(other_target_table['phot_rp_mean_mag'])
print("other G error ")
print(other_target_table['phot_g_mean_flux_error'])
print('other BP error')
print(other_target_table['phot_bp_mean_flux_error'])
print('other RP error')
print(other_target_table['phot_rp_mean_flux_error'])
print('other parallax')
print(other_target_table['parallax'])
#print('other parallax filled')
#print(other_target_table['parallax'].filled())
print('other parallax error')
print(other_target_table['parallax_error'])
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
    #rp_dist= remove_negative(rp_dist)
    #bp_dist= remove_negative(bp_dist)
    bp_mag_dist = get_mag(bp_dist, 'bp')
    rp_mag_dist = get_mag(rp_dist, 'rp')
    #make them be the same size
    smaller= np.min([bp_mag_dist.shape[0], rp_mag_dist.shape[0]])
    bp_mag_dist= bp_mag_dist[:smaller]
    rp_mag_dist= rp_mag_dist[:smaller]
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




#def correct_bp_flux(table, teff= teff, logg=logg):
    #obs_bp_flux, obs_bp_flux_dist = get_filter_vals(table, 'bp')
    #obs_rp_flux, obs_rp_flux_dist = get_filter_vals(table, 'rp')
    
    #return

q_dist = get_mc_distribution(q, q_err)

def calc_ns_mass(comp_mass, q=q):
    try:
        comp_mass, q = match_sizes(comp_mass, q)
    except IndexError:
        pass
    return q*comp_mass

def calc_mean_density(mass, radius):
    return (mass/(4./3 *np.pi*radius**3)).cgs
#def get_mass(radius, logg=logg):
    #logg = logg*u.
    
def get_comp_mass(q=q, m_psr = 1.4*u.Msun):
    """
    Get the companion mass by using the mass ratio and a pulsar mass assumption
    """
    return m_psr/q




def calc_min_period(mass, radius):
    """
    Returns the minimum period in hours that a given density star could have. 
    immediately converted to days though, so it actually returns the value in days
    """
    mean_density = calc_mean_density(mass, radius).cgs.value
    return ((107./mean_density)**2*u.hour).to(u.day)

def get_all_extinctions(table, logg= logg, teff=teff):
    """
    Return all each of the A_G values for the target using each of the methods, so they can all be compared.
    """
    gmag = table['phot_g_mean_mag']
    rpmag= table['phot_rp_mean_mag']
    model_g_absmag= pmc.get_model_absmag(logg = logg, teff= teff, passband_string = 'g')
    model_rp_absmag= pmc.get_model_absmag(logg= logg, teff = teff, passband_string = 'rp')
    koester_a_g=gaia_extinction.get_koester_a_g(model_g_absmag, model_rp_absmag, gmag, rpmag)
    obs_bp_rp = table['bp_rp']
    model_bp_rp = pmc.get_model_bp_rp(logg=logg, teff= teff)
    ebv_a_g = gaia_extinction.get_a_x(obs_bp_rp, model_bp_rp, passband_string = 'g')
    stat_a_g = gaia_extinction.get_statistical_a_g(obs_bp_rp, model_bp_rp)
    return np.hstack([koester_a_g, ebv_a_g, stat_a_g])
    
    

def get_pass_abs_mag(table, plot_all = False, passband_string= 'g', verbose = True, use_extinction= False, logg=logg, teff=teff, ext_method = 'B-V', get_extinction = False):
    mean_flux, flux_dist = get_filter_vals(table, passband_string)
    mag = get_mag(mean_flux, passband_string)
    #mag_dist = get_mag(mean_flux, passband_string )
    flux_dist= remove_negative(flux_dist, verbose=verbose)
    mag_dist= get_mag(flux_dist, passband_string)
    if use_extinction:
        if ext_method== 'B-V':
            obs_bp_rp = table['bp_rp']
            model_bp_rp = pmc.get_model_bp_rp(logg=logg, teff= teff)
            a_x = gaia_extinction.get_a_x(obs_bp_rp, model_bp_rp, passband_string = passband_string)
            print("Extinction ", "a_"+passband_string+":", a_x)
        elif ext_method== 'MG-MRP':
            #if (passband_string != 'g'):
                #print "*************************"
                #print "WARNING!!! YOU'RE USING THE G EXTINCTION FOR A DIFFERENT PASSBAND: "  +passband_string + " !!!!"
                #print "***************************"
            gmag = table['phot_g_mean_mag']
            rpmag= table['phot_rp_mean_mag']
            model_g_absmag= pmc.get_model_absmag(logg = logg, teff= teff, passband_string = 'g')
            model_rp_absmag= pmc.get_model_absmag(logg= logg, teff = teff, passband_string = 'rp')
            a_x=gaia_extinction.get_koester_a_g(model_g_absmag, model_rp_absmag, gmag, rpmag)
            print("Extinction ", "a_"+passband_string+":", a_x[0])
        elif ext_method == 'statistical':
            obs_bp_rp = table['bp_rp']
            model_bp_rp = pmc.get_model_bp_rp(logg=logg, teff= teff)
            a_x = gaia_extinction.get_statistical_a_g(obs_bp_rp, model_bp_rp)
        mag = mag-a_x
        mag_dist = mag_dist-a_x
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
    #distance_err= get_errors(distance_dist)
    #print("Distance:", distance[0], "-/+", distance_err[0,0], distance_err[1,0])
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
    if get_extinction:
        return abs_mag, abs_mag_error, abs_mag_dist, a_x
    else:
        return abs_mag, abs_mag_error, abs_mag_dist

def get_rad_mass(table, teff=teff, logg=logg, passband_string= 'g', plot_all = False, verbose= True, use_extinction= False, photometric= True, get_extinction = False, ext_method = 'B-V'):
    """
    Get the radius and mass for a given model when provided with the given absolute magnitude
    for a given band. Return the radius and mass for that absolute magnitude (or distribution) for the
    provided model.
    
    If photometric, then the radius is computed photometrically
    """
    #abs_mag,abs_mag_err, abs_mag_dist = get_pass_abs_mag(table, plot_all = plot_all, passband_string=passband_string, verbose = verbose)
    if photometric:
        if get_extinction:
            abs_mag,abs_mag_err, abs_mag_dist, a_x = get_pass_abs_mag(table, plot_all = plot_all, passband_string=passband_string, verbose = verbose, use_extinction= use_extinction, teff= teff, logg=logg, get_extinction= get_extinction, ext_method= ext_method)
        else:
            abs_mag,abs_mag_err, abs_mag_dist = get_pass_abs_mag(table, plot_all = plot_all, passband_string=passband_string, verbose = verbose, use_extinction= use_extinction, teff= teff, logg=logg, get_extinction= get_extinction, ext_method= ext_method)
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
    else:
        if get_extinction:
            abs_mag,abs_mag_err, abs_mag_dist, a_x = get_pass_abs_mag(table, plot_all = plot_all, passband_string=passband_string, verbose = verbose, use_extinction= use_extinction, teff= teff, logg=logg, ext_method= 'MG-MRP', get_extinction= get_extinction)
        else:
            abs_mag,abs_mag_err, abs_mag_dist = get_pass_abs_mag(table, plot_all = plot_all, passband_string=passband_string, verbose = verbose, use_extinction= use_extinction, teff= teff, logg=logg, ext_method= 'MG-MRP', get_extinction= get_extinction)
        #radius = pmc.get_radius(abs_mag, teff= teff, logg= logg, passband_string= passband_string)
        #mass = (pmc.get_mass(radius, logg)).to(u.Msun)
        #mass = get_comp_mass()
        #mass_dist= get_comp_mass(q=q_dist)
        mass = get_comp_mass(m_psr= m_psr)
        mass_dist= get_comp_mass(q=q_dist, m_psr= m_psr)
        radius = pmc.get_radius_from_mass(mass, logg)
        radius_dist= pmc.get_radius_from_mass(mass_dist, logg)
        #radius_dist = pmc.get_radius(abs_mag_dist, teff = teff, logg=logg, passband_string = passband_string)
        #mass_dist = (pmc.get_mass(radius_dist, logg)).to(u.Msun)
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
    ns_mass_dist= calc_ns_mass(mass_dist,q=q_dist) #now that we have a distribution for the mass ratio to be as well
    ns_mass_err= get_errors(ns_mass_dist)
    min_period = calc_min_period(mass, radius)
    min_period_dist= calc_min_period(mass_dist, radius_dist)
    min_period_err = get_errors(min_period_dist)
    print(passband_string+ " Abs Mag:", abs_mag[0], "-/+", abs_mag_err[0,0], abs_mag_err[1,0])
    try:
        print(passband_string+ " Radius:", radius[0], "-/+", radius_err[0,0], radius_err[1,0])
        print(passband_string+ " Comp Mass:", mass[0], "-/+", mass_err[0,0],mass_err[1,0])
        print(passband_string+ " PSR Mass:",ns_mass[0], "-/+", ns_mass_err[0,0],ns_mass_err[1,0])
        print(passband_string+ " Mean Density:", mean_density[0], "-/+", mean_density_err[0,0], mean_density_err[1,0])
        print(passband_string+ " Min Period:", min_period[0], "-/+", min_period_err[0,0], min_period_err[1,0])
    except TypeError:
        print(passband_string+ " Radius:", radius, "-/+", radius_err[0,0], radius_err[1,0])
        print(passband_string+ " Comp Mass:", mass, "-/+", mass_err[0,0],mass_err[1,0])
        print(passband_string+ " PSR Mass:",ns_mass, "-/+", ns_mass_err[0,0],ns_mass_err[1,0])
        print(passband_string+ " Mean Density:", mean_density, "-/+", mean_density_err[0,0], mean_density_err[1,0])
        print(passband_string+ " Min Period:", min_period, "-/+", min_period_err[0,0], min_period_err[1,0])
    if get_extinction:
        output_element= np.hstack([teff, logg, abs_mag, abs_mag_err.T[0], radius, radius_err.T[0], mass, mass_err.T[0], ns_mass, ns_mass_err.T[0], mean_density, mean_density_err.T[0], min_period, min_period_err.T[0],a_x]).value
        output_header = start_header + ', A_'+passband_string
    else:
        output_element= np.hstack([teff, logg, abs_mag, abs_mag_err.T[0], radius, radius_err.T[0], mass, mass_err.T[0], ns_mass, ns_mass_err.T[0], mean_density, mean_density_err.T[0], min_period, min_period_err.T[0]]).value
        output_header= start_header
    #return radius, mass, radius_dist, mass_dist
    #return radius, mass, radius_dist, mass_dist, output_element
    return radius, mass, radius_dist, mass_dist, output_element, output_header



def plot_target_lists(target_lists):
    for target_list in target_lists:
        print(target_list)
        list_table = Table.read(target_list)
        list_table= list_table.filled()
        counter=0
        for row in list_table:
            print('===========')
            #print(row)
            other_target_g_absmag, other_target_g_absmag_err, other_target_g_absmag_dist= get_pass_abs_mag(row, plot_all = False)
            other_target_bp_rp, other_target_bp_rp_err = get_bp_rp(row, plot_all = False)
            #blue_ones= np.where(other_target_bp_rp<0)
            #dim_ones= np.where(other_target_g_absmag>11)
            #print(
            if ((other_target_g_absmag > 11) and ( other_target_bp_rp<0)):
                print("\n\n************* This one***********\n\n")
            if row['phot_bp_mean_mag']>1e18:
                other_target_bp_rp= bp_rp_fill
                other_target_bp_rp_err= 0
            if row['parallax']>1e18:
                other_target_g_absmag=g_abs_mag_fill
            else:
                pass
            if counter==0:
                if error_bar:
                    plt.errorbar(np.copy(other_target_bp_rp), np.copy(other_target_g_absmag), yerr = np.copy(other_target_g_absmag_err),  xerr = np.copy(other_target_bp_rp_err), marker = 'o', markersize = 6, color = pod.plot_dict[target_list]['color'], capsize = 4, label = pod.plot_dict[target_list]['label'], linestyle ='none')
                else:
                    plt.plot(np.copy(other_target_bp_rp), np.copy(other_target_g_absmag), marker = 'o', markersize = 8, color = pod.plot_dict[target_list]['color'], label = pod.plot_dict[target_list]['label'], linestyle ='none')
            else:
                if error_bar:
                    plt.errorbar(np.copy(other_target_bp_rp), np.copy(other_target_g_absmag), yerr = np.copy(other_target_g_absmag_err),  xerr = np.copy(other_target_bp_rp_err), marker = 'o', markersize = 6, color = pod.plot_dict[target_list]['color'], capsize = 4, linestyle ='none')
                else:
                    plt.plot(np.copy(other_target_bp_rp), np.copy(other_target_g_absmag), marker = 'o', markersize = 8, color = pod.plot_dict[target_list]['color'], linestyle ='none')
            print(other_target_bp_rp, other_target_g_absmag)
            print(row['name'])
            #plt.annotate(str(row['name']),xy=(np.copy(other_target_bp_rp), np.copy(other_target_g_absmag)), xycoords='data', xytext=(np.copy(other_target_bp_rp+0.05),np.copy(other_target_g_absmag-0.1)), textcoords= 'data' , fontsize=8, color =pod.plot_dict[target_list]['color'])
            counter+=1
    
    
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




#target_g_absmag, target_g_absmag_err, target_g_absmag_dist = get_pass_abs_mag(target_table, plot_all = True, passband_string = 'g', use_extinction=True)
#target_bp_rp, target_bp_rp_err= get_bp_rp(target_table, plot_all = True)

target_g_absmag, target_g_absmag_err, target_g_absmag_dist = get_pass_abs_mag(target_table, plot_all = True, passband_string = 'g', use_extinction=False)
ptarget_g_absmag, ptarget_g_absmag_err, ptarget_g_absmag_dist = get_pass_abs_mag(target_table, plot_all = True, passband_string = 'g', use_extinction=True)
print("raw Abs Mag:", target_g_absmag[0], "-/+", target_g_absmag_err[0,0], target_g_absmag_err[1,0])
print("ext corr Abs Mag:", ptarget_g_absmag[0], "-/+", ptarget_g_absmag_err[0,0], ptarget_g_absmag_err[1,0])
target_bp_rp, target_bp_rp_err= get_bp_rp(target_table, plot_all = True)




#target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all = True)
#target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all = True, use_extinction=True)
#target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist, tossaway_element = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all = True, use_extinction=True)
target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist, tossaway_element, nothing_header = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all = True, use_extinction=True)


target_g_radius_err = get_errors(target_g_radius_dist)
target_g_mass_err = get_errors(target_g_mass_dist)


sim_target_gabsmag, sim_target_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff= teff, radius =target_g_radius)
sim_target_gabsmag_dist, sim_target_bp_rp= pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = target_g_radius_dist)
sim_target_gabsmag_err = get_errors(sim_target_gabsmag_dist)


#sim_target_gabsmag, sim_target_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff= teff, radius =target_g_radius)
#sim_target_gabsmag_dist, sim_target_bp_rp= pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = target_g_radius_dist)
#sim_target_gabsmag_err = get_errors(sim_target_gabsmag_dist)


print("Target Radius:", target_g_radius, "-/+", target_g_radius_err)
print("Target Mass:", target_g_mass, "-/+", target_g_mass_err)



#model_g_absmag, model_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = test_radii)


########3 other

###### BP

#target_bp_radius, target_bp_mass,target_bp_radius_dist, target_bp_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'bp', plot_all = True)
#target_bp_radius, target_bp_mass,target_bp_radius_dist, target_bp_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'bp', plot_all = True, use_extinction= True)
target_bp_radius, target_bp_mass,target_bp_radius_dist, target_bp_mass_dist, tossaway_element, nothing_header = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'bp', plot_all = True, use_extinction= True)


#target_bp_radius_dist, target_bp_mass_dist = get_rad_mass(, logg=logg, teff= teff, passband_string= 'bp')
target_bp_radius_err = get_errors(target_bp_radius_dist)
target_bp_mass_err = get_errors(target_bp_mass_dist)

print("Target Radius bp:", target_bp_radius, "-/+", target_bp_radius_err)
print("Target Mass bp:", target_bp_mass, "-/+", target_bp_mass_err)


sim_target_bp_absmag, sim_target_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff= teff, radius =target_bp_radius)
sim_target_bp_absmag_dist, sim_target_bp_rp= pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = target_bp_radius_dist)
sim_target_bp_absmag_err = get_errors(sim_target_bp_absmag_dist)

#####3 RP



#target_rp_radius, target_rp_mass,target_rp_radius_dist, target_rp_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'rp', plot_all = True)
#target_rp_radius, target_rp_mass,target_rp_radius_dist, target_rp_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'rp', plot_all = True, use_extinction=True)
target_rp_radius, target_rp_mass,target_rp_radius_dist, target_rp_mass_dist, tossaway_element, nothing_header = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'rp', plot_all = True, use_extinction=True)

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


#the main point
#uncomment this to actually plot the pulsar 2018-12-30
#plt.errorbar(target_bp_rp, target_g_absmag, yerr = target_g_absmag_err, xerr = target_bp_rp_err, marker = 'o', markersize = 12, color = 'g', capsize = 4, label = target_label, linestyle = 'none')

if single_list:
    for row in other_target_table:
        print('===========')
        #print(row)
        other_target_g_absmag, other_target_g_absmag_err, other_target_g_absmag_dist= get_pass_abs_mag(row, plot_all = False)
        other_target_bp_rp, other_target_bp_rp_err = get_bp_rp(row, plot_all = False)
        #other_target_bp_rp.fill_value=bp_rp_fill
        #other_target_g_absmag.fill_value= g_abs_mag_fill
        #other_target_bp_rp.filled()
        #other_target_g_absmag.filled()
        #print('other bp_rp',other_target_bp_rp)
        #print(other_target_bp_rp >"--")
        #print(other_target_bp_rp == np.nan)
        #if  other_target_bp_rp=='--':
            #print('it was')
            #other_target_bp_rp= 1.5
            #other_target_bp_rp_err= 1.5
        if row['phot_bp_mean_mag']>1e18:
            other_target_bp_rp= bp_rp_fill
            other_target_bp_rp_err= 0
        if row['parallax']>1e18:
            other_target_g_absmag=g_abs_mag_fill
        plt.errorbar(np.copy(other_target_bp_rp), np.copy(other_target_g_absmag), yerr = np.copy(other_target_g_absmag_err),  xerr = np.copy(other_target_bp_rp_err), marker = 'o', markersize = 8, color = list_color, capsize = 4, label = other_target_label, linestyle ='none')
        print(other_target_bp_rp, other_target_g_absmag)
        print(row['name'])
        #plt.text(np.copy(other_target_bp_rp+0.1), np.copy(other_target_g_absmag-0.05), str(row['dec']), fontsize=8, color ='b')
        ##plt.annotate(str(row['dec']),xy=(np.copy(other_target_bp_rp), np.copy(other_target_g_absmag)), xycoords='data', xytext=(np.copy(other_target_bp_rp+0.1),np.copy(other_target_g_absmag-0.1)), textcoords= 'data',arrowprops=dict(arrowstyle="<->") , fontsize=8, color ='b')
        plt.annotate(str(row['name']),xy=(np.copy(other_target_bp_rp), np.copy(other_target_g_absmag)), xycoords='data', xytext=(np.copy(other_target_bp_rp+0.1),np.copy(other_target_g_absmag-0.1)), textcoords= 'data' , fontsize=8, color =list_color)
else:
    plot_target_lists(pod.target_lists)
    plt.legend()
#plt.errorbar(other_target_bp_rp, other_target_g_absmag, yerr = other_target_g_absmag_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')
####This one \/ \/ \/
####plt.errorbar(other_target_bp_rp, other_target_g_absmag, yerr = other_target_g_absmag_err, xerr= other_target_bp_rp_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')

####plt.errorbar(other_target_bp_rp, other_target_g_absmag, xerr= other_target_bp_rp_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')

####plt.plot(model_bp_rp*np.ones(model_g_absmag.shape[0]), model_g_absmag, marker = '*', markersize= 8, color = 'green', label = "Model spectra at various radii " + str(np.round(test_radii.value, precision))+ " logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')

#uncomment this to plot the corrected version of the pulsar (companion) 2018-12-30
#plt.errorbar(sim_target_bp_rp, sim_target_gabsmag, yerr = sim_target_gabsmag_err, xerr = target_bp_rp_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = "Corrected", linestyle = 'none')

#plt.plot(sim_target_bp_rp, sim_target_gabsmag, marker = '*', markersize= 8, color = 'green', label = "R_G= "+str(np.round(target_g_radius, precision)[0]) + " M= "+ str(np.round(target_g_mass, precision)[0])+" logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')

#plt.plot(sim_target_bp_rp, sim_target_bp_absmag, marker = '*', markersize= 8, color = 'cyan', label = "R_BP= "+str(np.round(target_bp_radius, precision)[0]) + " M= "+ str(np.round(target_bp_mass, precision)[0])+" logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')

#plt.plot(sim_target_bp_rp, sim_target_rp_absmag, marker = '*', markersize= 8, color = 'magenta', label = "R_RP= "+str(np.round(target_rp_radius, precision)[0]) + " M= "+ str(np.round(target_rp_mass, precision)[0])+" logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')

###plt.plot(model_bp_rp*np.ones(model_g_absmag.shape[0]), model_g_absmag, marker = '*', markersize= 8, color = 'cyan', label = "Model spectra at various radii  logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')
###plt.legend()


#plt.plot( mass_bp_rp*np.ones(mass_gabsmag.shape[0]), mass_gabsmag,marker = '*', markersize= 8, color = 'magenta', label = "using the mass ratio and surface gravity" + str(logg) + " Teff: " +str(teff), linestyle = 'none')
##plt.legend()
#make_density_plot( generic_g_absmag,generic_bp_rp)
#plt.gca().invert_yaxis()

#plt.show()

#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
counts = polything.get_array()
print(counts.shape)
counts= np.sqrt(counts)
polything.set_array(counts)
polything.autoscale()

#plt.title('PSR J1431-4715' + title_suffix)

#plt.ylim([-4, 16])
plt.ylim(axes_y)

plt.gca().invert_yaxis()
plt.xlabel(r'$G_{BP} - G_{RP}$')
#plt.xlim([-1,5])
plt.xlim(axes_x)
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
teffs_written = []
teffs_wanted = np.arange(good_teff_range[0], good_teff_range[1]+1000, 1000)
loggs_written = []
loggs_wanted= np.arange(3, 7, 1)
output_list= []
#model_Rc=[]
#model_Rc_lo=[]
#model_Rc_hi=[]
#model_Mc_lo=[]
#model_Mc_hi
#model_Mpsr=[]
#model_Mpsr_lo=[]

#model_
for logg,teff in zip( logg_array,teff_array):

    print("Teff:", teff, "logg:", logg)
    model = wd(Teff = teff , logg = logg)
    #target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all= False, verbose= False)
    #target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all= False, verbose= False, use_extinction= True)
    #target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all= False, verbose= False, use_extinction= True, photometric= False)
    #target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist, output_element, output_header = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all= False, verbose= False, use_extinction= True, photometric= False, get_extinction=True)
    #target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist, output_element, output_header = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all= False, verbose= False, use_extinction= True, photometric= True, get_extinction=True)
    target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist, output_element, output_header = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all= False, verbose= False, use_extinction= True, photometric= True, get_extinction=True, ext_method = ext_method)
    #target_g_radius, target_g_mass,target_g_radius_dist, target_g_mass_dist = get_rad_mass(target_table, logg=logg, teff= teff, passband_string= 'g', plot_all= False, verbose= False, use_extinction= False)

    target_g_radius_err = get_errors(target_g_radius_dist)
    target_g_mass_err = get_errors(target_g_mass_dist)
    sim_target_gabsmag, sim_target_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff= teff, radius =target_g_radius)
    sim_target_gabsmag_dist, sim_target_bp_rp= pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = target_g_radius_dist)
    sim_target_gabsmag_err = get_errors(sim_target_gabsmag_dist)
    try:
        print("Model g Abs Mag:", sim_target_gabsmag[0], "-/+", sim_target_gabsmag_err[0,0], sim_target_gabsmag_err[1,0])
    except IndexError:
        print("Model g Abs Mag:", sim_target_gabsmag, "-/+", sim_target_gabsmag_err[0,0], sim_target_gabsmag_err[1,0])
    new_element= np.hstack([sim_target_gabsmag, sim_target_gabsmag_err.T[0], sim_target_bp_rp])
    extinction_values = get_all_extinctions(target_table, logg=logg, teff=teff)
    #total_element= np.hstack([output_element, new_element])
    total_element= np.hstack([output_element, new_element, extinction_values])
    output_list.append(total_element)
    print( "BP-RP:", sim_target_bp_rp)
    #print("G absmag:", sim_target_gabsmag, "-/+", sim_target_gabsmag_err[0,0], sim_target_gabsmag_err[1,0])
    #print("comp_mass:", target_g_mass, "-/+", get_errors(target_g_mass_dist))
    #print("NS_mass: ", calc_ns_mass(target_g_mass), "-/+", get_errors(calc_ns_mass(target_g_mass_dist)))
    print("--------------------------")
    #if ((teff.value >= good_teff_range[0]) and (teff.value <= good_teff_range[1])):
        #plt.errorbar(sim_target_bp_rp, sim_target_gabsmag, yerr= sim_target_gabsmag_err, marker = '*', markersize= 8, color = 'green')
        ##if ((int(teff.value) not in teffs_written) and (int(teff.value) in teffs_wanted)):
        #if ((int(teff.value) not in teffs_written) and (int(teff.value) in teffs_wanted) and (logg == 5.0)):
            ##plt.text(sim_target_bp_rp, sim_target_gabsmag, str(teff.value)+"K", rotation = 45, horizontalalignment = 'center', verticalalignment= 'center')
            #if teff.value == 10000:
                #plt.annotate( str(teff.value)+"K", xy= (sim_target_bp_rp, sim_target_gabsmag), xycoords = 'data', xytext = (sim_target_bp_rp - 0.5, sim_target_gabsmag +1), rotation = 45, arrowprops=dict(facecolor='black', shrink=0.05))
            #else:
                #plt.annotate( str(teff.value)+"K", xy= (sim_target_bp_rp, sim_target_gabsmag), xycoords = 'data', xytext = (sim_target_bp_rp - 0.5, sim_target_gabsmag +2), rotation = 45, arrowprops=dict(facecolor='black', shrink=0.05))
            #teffs_written.append(int(teff.value))
    if ((teff.value >= good_teff_range[0]) and (teff.value <= good_teff_range[1])):
        plt.errorbar(sim_target_bp_rp, sim_target_gabsmag, yerr= sim_target_gabsmag_err, marker = '*', markersize= 8, color = 'green')
        #if ((int(teff.value) not in teffs_written) and (int(teff.value) in teffs_wanted)):
        if ((logg not in loggs_written) and (logg in loggs_wanted) and (teff.value== teffs_wanted.min())):
            #plt.text(sim_target_bp_rp, sim_target_gabsmag, str(teff.value)+"K", rotation = 45, horizontalalignment = 'center', verticalalignment= 'center')
            #if teff.value == 10000:
                #plt.annotate( str(teff.value)+"K", xy= (sim_target_bp_rp, sim_target_gabsmag), xycoords = 'data', xytext = (sim_target_bp_rp - 0.5, sim_target_gabsmag +1), rotation = 45, arrowprops=dict(facecolor='black', shrink=0.05))
            #plt.annotate( "log(g)="+ str(logg), xy= (sim_target_bp_rp, sim_target_gabsmag), xycoords = 'data', xytext = (sim_target_bp_rp - 0.6, sim_target_gabsmag), arrowprops=dict(facecolor='black', shrink=0.05))
            plt.annotate( "log(g)="+ str(logg), xy= (sim_target_bp_rp, sim_target_gabsmag), xycoords = 'data', xytext = (sim_target_bp_rp +2, sim_target_gabsmag-2), arrowprops=dict(facecolor='black', shrink=0.05))
            loggs_written.append(logg)

output_header= output_header+ ", model_g_absmag, model_g_absmag_err_lo, model_g_absmag_err_hi, model_bp_rp"

output_header = output_header + ', a_g_koester, a_g_e_b_v, a_g_stat'
output_array = np.array(output_list)
np.savetxt(output_filename, output_array, delimiter= ',', header= output_header)
plt.errorbar(target_bp_rp, target_g_absmag, yerr = target_g_absmag_err, xerr = target_bp_rp_err, marker = '*', markersize = 12, color = 'b', capsize = 4, label = target_label, linestyle = 'none')










#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
counts = polything.get_array()
print(counts.shape)
counts= np.sqrt(counts)
polything.set_array(counts)
polything.autoscale()

plt.title('PSR J1431-4715' + title_suffix)

#plt.ylim([-4, 16])
plt.ylim(axes_y)

plt.gca().invert_yaxis()
plt.xlabel(r'$G_{BP} - G_{RP}$')
#plt.xlim([-1,5])
plt.xlim(axes_x)
plt.ylabel(r'$M_G$')
#plt.legend()
plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)
plt.show()






