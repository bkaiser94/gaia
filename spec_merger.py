"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-03-31

Take two SDSS spectra that also have Gaia parallaxes and stuff, and then rescale the flux in the spectra to then be at the same distance and add the fluxes together to produce a hybrid spectrum.

Should also plot the rescaled spectrum of each individual star at the same time as a quality check and to demonstrate the relative contributions




"""
from __future__ import print_function
import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, QTable
import matplotlib.pyplot as plt
import scipy.stats as scistats
from astropy.io import fits as iofits
#import seaborn as sns
import astropy

import gaia_extinction
#import wdatmos
import plotting_dicts as pod
import spec_plot_tools as spt
########################################
########################################
########################################
########################################
############

parallax_correction = 0.029 #from Lindgren et al 2018



#mdwarf_spec_file= 'sdssj1246p3608_sdss_spec.fits'
#mdwarf_gaia_file='sdssj1246p3608_gaia.csv'

#mdwarf_spec_file='sdssj1408p2021_sdss_spec.fits'
#mdwarf_gaia_file='sdssj1408p2021_gaia.csv'


#mdwarf_spec_file='2MASSJ1055p0808_sdss_spec.fits'
#mdwarf_gaia_file='2MASSJ1055p0808_gaia.csv'


mdwarf_spec_file='2MASSJ1458p2839_sdss_spec.fits'
mdwarf_gaia_file='2MASSJ1458p2839_gaia.csv'

wdwarf_spec_file='WD1401p457_sdss_spec.fits'
wdwarf_gaia_file='20190218_test_ultcool.csv'

wd_name= wdwarf_spec_file.split('_')[0]
mdwarf_name= mdwarf_spec_file.split('_')[0]

wdwarf_gaia_table= Table.read(wdwarf_gaia_file)[0]
mdwarf_gaia_table= Table.read(mdwarf_gaia_file)[0]

wdwarf_spec_hdu= iofits.open(wdwarf_spec_file)
mdwarf_spec_hdu= iofits.open(mdwarf_spec_file)

wdwarf_spec_array= wdwarf_spec_hdu[1].data
print(wdwarf_spec_array)
print(wdwarf_spec_array['loglam'])

mdwarf_spec_array = mdwarf_spec_hdu[1].data

#print(wdwarf_spec_array['loglam']-mdwarf_spec_array['loglam'])


###################################
#################################
#################################
#################################


def plot_raw_spec(spec_hdu, name=''):
    plt.plot(spec_hdu['loglam'], spec_hdu['flux'])
    plt.title(name)
    plt.show()
    return



def distance_modulus(g_mag, distance):
    return g_mag - 5*np.log10(distance/10.)
    

def get_mag(flux):
    return -2.5*np.log10(flux) 

def mag_to_flux(mag):
    return 10.**((mag)/-2.5)


def get_abs_flux(table, spec_hdu,  plot_all = False,  verbose = True):
    #mean_flux, flux_dist = get_filter_vals(table, passband_string)
    mag = get_mag(spec_hdu['flux'])
    #flux_dist= remove_negative(flux_dist, verbose=verbose)
    #mag_dist= get_mag(flux_dist, passband_string)
    parallax = table['parallax']+parallax_correction
    parallax = parallax*1e-3
    distance = 1./parallax
    parallax_error = table['parallax_error']*1e-3
    #parallax_dist = get_mc_distribution(parallax, parallax_error)
    #parallax_dist = remove_negative(parallax_dist, verbose= verbose)
    #if verbose:
        #print(passband_string+ "_calc" + "-" + passband_string+ "_measured", mag - table['phot_' +passband_string+'_mean_mag'])
    #else:
        #pass
    #if parallax < 0:
        #parallax_median = np.nanmedian(parallax_dist)
        #if verbose:
            #print("PARALLAX < 0!", parallax, "setting to median of positive distribution:", parallax_median)
        #parallax = parallax_median
    #else:
        #pass
    distance = 1./parallax
    #distance_dist = 1./parallax_dist
    #index_length = distance_dist.shape[0]
    #print("index_length",index_length)
    #print("mag_dist.shape", mag_dist.shape)
    #mag_dist = mag_dist[:index_length]
    abs_mag = distance_modulus(mag, distance)
    #abs_mag_dist = distance_modulus(mag_dist, distance_dist)
    #abs_mag_error= get_errors(abs_mag_dist)
    #if plot_all:
        #plt.hist(abs_mag_dist, bins=75, normed=1, label = 'MC Distribution', color = 'g')
        #plt.axvline(np.nanmedian(abs_mag_dist), color = 'k', linestyle = '--', label = 'Median of MC Dist')
        #plt.axvline(np.nanpercentile(abs_mag_dist, 84), color = 'cyan')
        #plt.errorbar(abs_mag, 0.5, xerr = abs_mag_error, marker = '*', markersize = 8, color = 'b', label = "M_"+passband_string, capsize = 4)
        #plt.xlabel('M_'+ passband_string)
        #plt.legend()
        #plt.show()
    #else:
        #pass
    abs_flux= mag_to_flux(abs_mag)
    return abs_flux

#################################
#################################
#################################
#################################
#################################



plot_raw_spec(wdwarf_spec_array, name=wd_name)
plot_raw_spec(mdwarf_spec_array, name=mdwarf_name)

### convolve the spectra to dampen noise
mdwarf_spec_array= spt.convolve_spec(mdwarf_spec_array, kernel_type= 'box', width = 2)
wdwarf_spec_array= spt.convolve_spec(wdwarf_spec_array, kernel_type= 'box', width = 2)


plot_raw_spec(wdwarf_spec_array, name=wd_name)
plot_raw_spec(mdwarf_spec_array, name=mdwarf_name)

wdwarf_spec_array= spt.interpolate_spec(mdwarf_spec_array, wdwarf_spec_array)


plot_raw_spec(wdwarf_spec_array, name=wd_name)
plot_raw_spec(mdwarf_spec_array, name=mdwarf_name)


wdwarf_abs_flux= get_abs_flux(wdwarf_gaia_table, wdwarf_spec_array)

mdwarf_abs_flux= get_abs_flux(mdwarf_gaia_table, mdwarf_spec_array)

plt.plot(10**wdwarf_spec_array['loglam'], wdwarf_abs_flux, label=wd_name)
plt.plot(10**mdwarf_spec_array['loglam'], mdwarf_abs_flux, label=mdwarf_name)
plt.legend()
plt.show()


#plt.plot(10**wdwarf_spec_array['loglam'], wdwarf_abs_flux, label='WD')
#plt.plot(10**mdwarf_spec_array['loglam'], mdwarf_abs_flux, label='dM')
#plt.legend()
#plt.show()

merged_flux= wdwarf_abs_flux+mdwarf_abs_flux


plt.plot(10**wdwarf_spec_array['loglam'], wdwarf_abs_flux, label=wd_name)
plt.plot(10**mdwarf_spec_array['loglam'], mdwarf_abs_flux, label=mdwarf_name)
plt.plot(10**wdwarf_spec_array['loglam'], merged_flux, label='merged')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel(r'Flux $[10^{-17}erg/cm^2/s/\AA]$')
plt.xlim(np.nanmin(10**wdwarf_spec_array['loglam']), np.nanmax(10**wdwarf_spec_array['loglam']))
plt.legend()
plt.show()


