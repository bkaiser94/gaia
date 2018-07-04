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
#import seaborn as sns
import astropy

import passband_model_convolution as pmc

plt.rc('font', size =18)
#plt.rc('lines', markersize=12)
#plt.rc('font', size = 11)
plt.rc('lines', markersize = 5)
precision = 2
parallax_correction = 0.029 #from Lindgren et al 2018
#parallax_correction = 0 #so nothing done

teff = 7250
logg = 6.0

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


zeropoint_dict={"g": [25.6884, 0.0018],
                "bp": [25.3514, 0.0014],
                "rp": [24.7619, 0.0019]} #from Evans et al 2018, the DR2 values [ZP, sigma]


#zeropoint_dict={"g": [ 25.7933969562,  0.0017848281],
                #"bp": [25.3805596387,  0.0013917453], 
                #"rp": [25.1161276701, 0.001914645] } #AB from Evans et al 2018, the DR2 values [ZP, sigma]
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
    try:
        return np.array([[median-low_bar],[high_bar-median]])
    except astropy.units.core.UnitsError as error:
        print(error)
        return  np.array([[median.value-low_bar],[high_bar-median.value]])
    #return np.vstack([np.array(median-low_bar),np.array(high_bar-median)])


def get_filter_vals(table, filter_string):
    phot_mean_flux = table['phot_'+filter_string+ '_mean_flux']
    phot_mean_flux_error = table['phot_' + filter_string + '_mean_flux_error']
    flux_distribution = get_mc_distribution(phot_mean_flux, phot_mean_flux_error)
    return phot_mean_flux, flux_distribution


def get_bp_rp(table, plot_all = False):
    bp_mean_flux, bp_dist = get_filter_vals(table, 'bp')
    rp_mean_flux, rp_dist = get_filter_vals(table, 'rp')
    bp_mag = get_mag(bp_mean_flux, 'bp')
    print("bp_calc-bp_measured", bp_mag - table['phot_bp_mean_mag'])
    rp_mag = get_mag(rp_mean_flux, 'rp')
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
    print("g_mag_dist.shape", g_mag_dist.shape, "distance_dist.shape", distance_dist.shape)
    g_mag_dist = g_mag_dist[:index_length]
    print("g_mag_dist.shape", g_mag_dist.shape, "distance_dist.shape", distance_dist.shape)
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


#def get_mass(radius, logg=logg):
    #logg = logg*u.
    
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

print(generic_g_absmag.shape)
#print(target_table['ra'])
#print(target_table['dec'])
#target_parallax = target_table['parallax']+parallax_correction
#target_parallax = target_parallax*1e-3
#target_distance = 1./target_parallax
#target_parallax_error = target_table['parallax_error']*1e-3
#target_parallax_dist = get_mc_distribution(target_parallax, target_parallax_error)
#target_parallax_dist = remove_negative(target_parallax_dist)
##remove_negative(target_parallax_dist)
#target_distance_dist = 1./target_parallax_dist
##target_distance =  1700 #from Jennings et al 2018 d_LK

##target_distance_err_bounds = 

#target_extinction= target_table['a_g_val']

#target_g_mag = target_table['phot_g_mean_mag']
##target_bp_rp = target_table['bp_rp']
#target_bp_rp, target_bp_rp_err = get_bp_rp(target_table, plot_all = True)


#target_g_flux = target_table['phot_g_mean_flux']


print("special thing: ", distance_modulus(18.44, 2.2*1000))



#generic_g_absmag = distance_modulus(generic_g_mag, generic_distance, extinction = generic_extinction)
#target_g_absmag = distance_modulus(target_g_mag, target_distance,extinction= target_extinction)
#target_g_absmag_dist = distance_modulus(target_g_mag, target_distance_dist, extinction = target_extinction)
#target_g_absmag_err = get_errors(target_g_absmag_dist)

#plt.hist(target_g_absmag_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
#plt.axvline(np.nanmedian(target_g_absmag_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
#plt.axvline(np.nanpercentile(target_g_absmag_dist, 84), color = 'cyan')
#plt.errorbar( target_g_absmag, 0.5, xerr = target_g_absmag_err, marker = '*', markersize = 8, color = 'b', label = target_label, capsize = 4)
##plt.axvline(x=target_g_absmag, color = 'r', linestyle = ':', label = 'Measured value')
#plt.xlabel(r'$M_G$')
#plt.title(target_label)
#plt.legend()
#plt.show()
target_g_absmag, target_g_absmag_err, target_g_absmag_dist = get_g_abs_mag(target_table, plot_all = True)
target_bp_rp, target_bp_rp_err= get_bp_rp(target_table, plot_all = True)

target_g_radius = pmc.get_radius(target_g_absmag, teff= teff, logg= logg, passband_string= 'G')
target_g_radius_dist = pmc.get_radius(target_g_absmag_dist, teff= teff, logg= logg, passband_string= 'G')


target_g_radius_err = get_errors(target_g_radius_dist)

target_g_mass = pmc.get_mass(target_g_radius, logg)
print("Mass: " , target_g_mass.to(u.Msun))

print("Radius for expected mass:", pmc.get_radius_from_mass(0.1325*u.Msun, logg))
print("Radius for double expected mass:", pmc.get_radius_from_mass(2*0.1325*u.Msun, logg))
print("Radius for triple expected mass:", pmc.get_radius_from_mass(3*0.1325*u.Msun, logg))
test_masses = 0.1325*np.array([1,2,3])*u.Msun
new_test_radii = pmc.get_radius_from_mass(test_masses, logg)
mass_gabsmag, mass_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff= teff, radius =new_test_radii)



plt.hist(target_g_radius_dist.value, bins=75, normed=1, label = 'MC Distribution', color = 'g')
plt.axvline(np.nanmedian(target_g_radius_dist.value), color = 'k', linestyle = '--', label = 'Median of MC Dist')
plt.axvline(np.nanpercentile(target_g_radius_dist.value, 84), color = 'cyan')
plt.errorbar(target_g_radius.value, 0.5, xerr = target_g_radius_err, marker = '*', markersize = 8, color = 'b', label = r"$R_{*}$", capsize = 4)
#plt.axvline(x=target_g_absmag, color = 'r', linestyle = ':', label = 'Measured value')
plt.xlabel(r'$R_{*}(R_{\odot}$')
#plt.title(target_label)
plt.legend()
plt.show()
print("Target Radius:", target_g_radius, "+/-", target_g_radius_err)


model_g_absmag, model_bp_rp = pmc.get_model_CMD_loc(logg= logg, teff = teff, radius = test_radii)


########3 other



#other_target_parallax = other_target_table['parallax']+parallax_correction
#other_target_parallax = other_target_parallax*1e-3
#other_target_distance = 1./other_target_parallax
####other_target_distance = 2.2*1000
###other_target_distance = 6.4*1000

#other_target_parallax_error = other_target_table['parallax_error']*1e-3
#other_target_parallax_dist = get_mc_distribution(other_target_parallax, other_target_parallax_error)
#other_target_parallax_dist = remove_negative(other_target_parallax_dist)
##remove_negative(target_parallax_dist)
#other_target_distance_dist = 1./other_target_parallax_dist

##target_distance_err_bounds = 

#other_target_extinction= other_target_table['a_g_val']

#other_target_g_mag =other_target_table['phot_g_mean_mag']
##other_target_bp_rp = other_target_table['bp_rp']
other_target_g_absmag, other_target_g_absmag_err, other_target_g_absmag_dist= get_g_abs_mag(other_target_table, plot_all = True)
other_target_bp_rp, other_target_bp_rp_err = get_bp_rp(other_target_table, plot_all = True)





#generic_g_absmag = distance_modulus(generic_g_mag, generic_distance, extinction = generic_extinction)
#other_target_g_absmag = distance_modulus(other_target_g_mag, other_target_distance)
#other_target_g_absmag_dist = distance_modulus(other_target_g_mag, other_target_distance_dist, extinction = other_target_extinction)
#other_target_g_absmag_err = get_errors(other_target_g_absmag_dist)

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

#plt.scatter(generic_bp_rp, generic_g_absmag, s = 1, alpha = 0.05, color = 'k')
#make_density_plot(generic_g_absmag, generic_bp_rp)
#plt.plot(target_bp_rp, target_g_absmag, marker = '*', markersize = 8, color = 'b')
#plt.errorbar(target_bp_rp, target_g_absmag, yerr = target_g_absmag_err, marker = '*', markersize = 8, color = 'b', capsize = 4, label = target_label, linestyle = 'none')
plt.errorbar(target_bp_rp, target_g_absmag, yerr = target_g_absmag_err, xerr = target_bp_rp_err, marker = '*', markersize = 8, color = 'b', capsize = 4, label = target_label, linestyle = 'none')

#plt.errorbar(other_target_bp_rp, other_target_g_absmag, yerr = other_target_g_absmag_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')
#This one \/ \/ \/
#plt.errorbar(other_target_bp_rp, other_target_g_absmag, yerr = other_target_g_absmag_err, xerr= other_target_bp_rp_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')

#plt.errorbar(other_target_bp_rp, other_target_g_absmag, xerr= other_target_bp_rp_err, marker = '*', markersize = 8, color = 'g', capsize = 4, label = other_target_label, linestyle ='none')

#plt.plot(model_bp_rp*np.ones(model_g_absmag.shape[0]), model_g_absmag, marker = '*', markersize= 8, color = 'green', label = "Model spectra at various radii " + str(np.round(test_radii.value, precision))+ " logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')

plt.plot(model_bp_rp*np.ones(model_g_absmag.shape[0]), model_g_absmag, marker = '*', markersize= 8, color = 'cyan', label = "Model spectra at various radii  logg: " + str(logg) + " Teff: " +str(teff), linestyle = 'none')
plt.legend()


plt.plot( mass_bp_rp*np.ones(mass_gabsmag.shape[0]), mass_gabsmag,marker = '*', markersize= 8, color = 'magenta', label = "using the mass ratio and surface gravity" + str(logg) + " Teff: " +str(teff), linestyle = 'none')
plt.legend()

polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
counts = polything.get_array()
print(counts.shape)
counts= np.sqrt(counts)
polything.set_array(counts)
polything.autoscale()
#plt.hist2d(generic_bp_rp, generic_g_absmag, bins = 100, cmap = 'Reds')
#plt.title('PSR J1431-4715' + title_suffix)
plt.title('PSR J1431-4715' + title_suffix)

plt.ylim([-4, 16])

plt.gca().invert_yaxis()
plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.xlim([-1,5])
plt.ylabel(r'$M_G$')
#plt.legend()
plt.subplots_adjust(wspace = 0, hspace = 0, top = 0.90, bottom = 0.10, left = 0.10, right = 0.90)
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


##################################3
def plot_names(name_array, g_absmag, bp_rp):
    for name, g, bprp  in zip(name_array, g_absmag, bp_rp):
        plt.text(bprp, g, name, fontsize = 8, color = 'b')

for thing in pulsar_list_all:
    print(thing)
pulsar_names = pulsar_list_all['PSR']
pulsar_names = np.genfromtxt('pulsar_names.txt', dtype = 'S32')
print(pulsar_names)
pulsar_g_mag= np.float_(pulsar_list_all['G_mag'])



##########################3

pulsar_abs_g_mag = distance_modulus(pulsar_g_mag, np.float_(pulsar_list_all['d_pi'])*1000)
pulsar_abs_g_mag_lo = distance_modulus(pulsar_g_mag, np.float_(pulsar_list_all['d_pi'])*1000-np.float_(pulsar_list_all['d_pi_lo'])*1000)
#pulsar_abs_g_mag_hi = distance_modulus(pulsar_g_mag, np.float_(pulsar_list_all['d_pi'])*1000+np.float_(pulsar_list_all['d_pi_hi'])*1000)
#pulsar_abs_g_mag_lo = pulsar_abs_g_mag-pulsar_abs_g_mag_lo #because higher numbers mean dimmer and farther, so 
#pulsar_abs_g_mag_hi = pulsar_abs_g_mag_hi- pulsar_abs_g_mag #difference on the high side
#pulsar_abs_g_mag_err = np.vstack([pulsar_abs_g_mag_lo, pulsar_abs_g_mag_hi]) #hopefully correctly shaped

standin_bp_rp = np.random.rand(pulsar_abs_g_mag.shape[0])
#bp_rp_errorbars = np.ones(standin_bp_rp.shape)*6
#plt.errorbar(standin_bp_rp, pulsar_abs_g_mag, yerr= pulsar_abs_g_mag_err, xerr= bp_rp_errorbars, linestyle = 'none', marker = '*', color = 'b', capsize = 4)
#plot_names(pulsar_names, pulsar_abs_g_mag, standin_bp_rp)
#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
#counts = polything.get_array()
#print(counts.shape)
#counts= np.sqrt(counts)
#polything.set_array(counts)
#polything.autoscale()
##plt.hist2d(generic_bp_rp, generic_g_absmag, bins = 100, cmap = 'Reds')
##plt.title('PSR J1431-4715' + title_suffix)
#plt.title('d_pi' + title_suffix +'(random BP-RP) from Jennings et al. 2018')

##plt.ylim([-4, 16])

#plt.gca().invert_yaxis()
#plt.xlabel(r'$G_{BP} - G_{RP}$')
##plt.xlim([-1,5])
#plt.ylabel(r'$M_G$')
##plt.legend()
#plt.show()



###############3



#pulsar_abs_g_mag = distance_modulus(pulsar_g_mag, np.float_(pulsar_list_all['d_LK'])*1000)
#pulsar_abs_g_mag_lo = distance_modulus(pulsar_g_mag, np.float_(pulsar_list_all['d_LK'])*1000-np.float_(pulsar_list_all['d_LK_lo'])*1000)
#pulsar_abs_g_mag_hi = distance_modulus(pulsar_g_mag, np.float_(pulsar_list_all['d_LK'])*1000+np.float_(pulsar_list_all['d_LK_hi'])*1000)
#pulsar_abs_g_mag_lo = pulsar_abs_g_mag-pulsar_abs_g_mag_lo #because higher numbers mean dimmer and farther, so 
#pulsar_abs_g_mag_hi = pulsar_abs_g_mag_hi- pulsar_abs_g_mag #difference on the high side
#pulsar_abs_g_mag_err = np.vstack([pulsar_abs_g_mag_lo, pulsar_abs_g_mag_hi]) #hopefully correctly shaped

##standin_bp_rp = np.random.rand(pulsar_abs_g_mag.shape[0])
#print(standin_bp_rp.shape)
#print(pulsar_abs_g_mag.shape)
##plt.errorbar(standin_bp_rp, pulsar_abs_g_mag, yerr= pulsar_abs_g_mag_err, xerr= bp_rp_errorbars, linestyle = 'none', marker = '*', color = 'b', capsize = 4)
#plt.errorbar(standin_bp_rp, pulsar_abs_g_mag, yerr= pulsar_abs_g_mag_err, linestyle = 'none', marker = '*', color = 'b', capsize = 4)
#plot_names(pulsar_names, pulsar_abs_g_mag, standin_bp_rp)

#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
#counts = polything.get_array()
#print(counts.shape)
#counts= np.sqrt(counts)
#polything.set_array(counts)
#polything.autoscale()
##plt.hist2d(generic_bp_rp, generic_g_absmag, bins = 100, cmap = 'Reds')
##plt.title('PSR J1431-4715' + title_suffix)
#plt.title('d_LK' + title_suffix + ' random BP-RP from Jennings et al 2018')

##plt.ylim([-4, 16])

#plt.gca().invert_yaxis()
#plt.xlabel(r'$G_{BP} - G_{RP}$')
##plt.xlim([-1,5])
#plt.ylabel(r'$M_G$')
##plt.legend()
#plt.show()


######################

#pulsar_parallax =  np.float_(pulsar_list_all['pi'])
#pulsar_parallax_error = np.float_(pulsar_list_all['pi_error'])

#pulsar_parallax =pulsar_parallax*1e-3
#pulsar_distance = 1./pulsar_parallax
#pulsar_parallax_error = pulsar_parallax_error*1e-3

#print(pulsar_parallax.shape, pulsar_parallax_error.shape)
#pulsar_parallax_dist = get_mc_distribution(pulsar_parallax, pulsar_parallax_error)
#pulsar_parallax_dist = remove_negative(pulsar_parallax_dist)
##remove_negative(target_parallax_dist)
#pulsar_distance_dist = 1./pulsar_parallax_dist

#pulsar_g_absmag = distance_modulus(pulsar_g_mag, pulsar_distance)
#pulsar_g_absmag_dist = distance_modulus(pulsar_g_mag, pulsar_distance_dist)
##pulsar_g_absmag_err = get_errors(pulsar_g_absmag_dist)
#print("pulsar_g_absmag_dist.shape", pulsar_g_absmag_dist.shape)
#pulsar_g_agmsag_err = np.std(pulsar_g_absmag_dist, axis = 0)


#standin_bp_rp = np.random.rand(pulsar_abs_g_mag.shape[0])
#plt.errorbar(standin_bp_rp, pulsar_g_absmag, yerr= pulsar_g_absmag_err, linestyle = 'none', marker = '*', color = 'b', capsize = 4)
#plot_names(pulsar_names, pulsar_g_absmag, standin_bp_rp)
#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(1000,1000), cmap = 'hot', mincnt = 1)
#polything = plt.hexbin(generic_bp_rp, generic_g_absmag, gridsize=(grid_num, grid_num), cmap = 'hot', mincnt = 1, label = "H-R")
#counts = polything.get_array()
#print(counts.shape)
#counts= np.sqrt(counts)
#polything.set_array(counts)
#polything.autoscale()
##plt.hist2d(generic_bp_rp, generic_g_absmag, bins = 100, cmap = 'Reds')
##plt.title('PSR J1431-4715' + title_suffix)
#plt.title('Inverse parallax' + title_suffix +'(random BP-RP)')

##plt.ylim([-4, 16])

#plt.gca().invert_yaxis()
#plt.xlabel(r'$G_{BP} - G_{RP}$')
##plt.xlim([-1,5])
#plt.ylabel(r'$M_G$')
##plt.legend()
#plt.show()






















