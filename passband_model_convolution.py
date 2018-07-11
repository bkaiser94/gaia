"""
Created by Ben Kaiser (UNC-Chapel Hill)

Plot the DR2 passbands and the "Revised DR2 passbands" for the GAIA mission.

Requires one to already have downloaded the gaia passband files from the gaia website and update the file paths to get to them.

"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from glob import glob
import scipy.optimize as sciop
from astropy.time import Time
from astropy import coordinates as coords
from astropy import units as u
from astropy import constants as const
from astropy import convolution as conv
import scipy.interpolate as scinterp



import wdatmos
import spec_plot_tools as spt

teff= 7250
logg = 6.0
wd=wdatmos.wdmodel(filename='ELM.hdf5')
model = wd(Teff = teff, logg = logg)

distance = 1.56*1000 #pc

radius = 0.1*u.Rsun

area = 0.725*u.m**2
#area= 1.

#model_waves = model['w'].data
#model_flux = model['flux'].data #since we'll be arbitrarily-ish scaling this it won't work.
#model_flux = model_flux

#model_spec  = np.vstack([model_waves, model_flux])

#This should really start coming from the actual text file that contains these values, and I should load them in here
zeropoint_dict={"g": [25.6884, 0.0018],
                "bp": [25.3514, 0.0014],
                "rp": [24.7619, 0.0019]} # Vega from Evans et al 2018, the DR2 values [ZP, sigma]


#zeropoint_dict={"G": [ 25.7933969562,  0.0017848281],
                #"BP": [25.3805596387,  0.0013917453], 
                #"RP": [25.1161276701, 0.001914645] } #AB from Evans et al 2018, the DR2 values [ZP, sigma]


dr2_passband_file = 'GaiaDR2_Passbands_ZeroPoints/GaiaDR2_Passbands.dat'
dr2_rev_passband_file = 'GaiaDR2_Revised_Passbands_ZeroPoints/GaiaDR2_RevisedPassbands.dat'


#dr2_all = np.genfromtxt(dr2_passband_file).T
dr2_all = np.genfromtxt(dr2_rev_passband_file).T #I have changed the passbands to use the revised transmission curves
#rev_all = np.genfromtxt(dr2_rev_passband_file).T
print(dr2_all.shape)
wavelengths_dr2 = dr2_all[0]*10 #converted to angstroms
#wavelengths_rev = rev_all[0]
print(wavelengths_dr2.shape)

def good_for_plots_dr2(bandpass, sigma):
    good_vals = np.where(bandpass < 99)
    return wavelengths_dr2[good_vals], bandpass[good_vals], sigma[good_vals]

#def good_for_plots_rev(bandpass, sigma):
    #good_vals = np.where(bandpass < 99)
    #return wavelengths_rev[good_vals], bandpass[good_vals], sigma[good_vals]


Gband_dr2 = dr2_all[1]
#Gband_rev = rev_all[1]

Gband_sig_dr2 = dr2_all[2]
#Gband_sig_rev = rev_all[2]

BPband_dr2 = dr2_all[3]
#BPband_rev = rev_all[3]

BPband_sig_dr2 = dr2_all[4]
#BPband_sig_rev = rev_all[4]

RPband_dr2 = dr2_all[5]
#RPband_rev = rev_all[5]

RPband_sig_dr2 = dr2_all[6]
#RPband_sig_rev = rev_all[6]

dr2_G_tuple = good_for_plots_dr2(Gband_dr2, Gband_sig_dr2)
#rev_G_tuple = good_for_plots_rev(Gband_rev, Gband_sig_rev)

dr2_RP_tuple = good_for_plots_dr2(RPband_dr2, RPband_sig_dr2)
#rev_RP_tuple = good_for_plots_rev(RPband_rev, RPband_sig_rev)

dr2_BP_tuple = good_for_plots_dr2(BPband_dr2, BPband_sig_dr2)
#rev_BP_tuple = good_for_plots_rev(BPband_rev, BPband_sig_rev)

#lower_wavebound = np.nanmin(dr2_BP_tuple[0])
#upper_wavebound = np.nanmax(dr2_RP_tuple[0])
#model_spec= spt.clean_spectrum(model_spec, lower_wavebound, upper_wavebound, [])


def get_model_spec(teff= teff, logg= logg):
    wd=wdatmos.wdmodel(filename='ELM.hdf5')
    model = wd(Teff = teff, logg = logg)
    model_waves = model['w'].data
    #model_flux = model['flux'].data #since we'll be arbitrarily-ish scaling this it won't work.
    #model_flux = model['flux'].data/4. #since we'll be arbitrarily-ish scaling this it won't work.
    model_flux =np.pi* model['flux'].data #since we'll be arbitrarily-ish scaling this it won't work.
    model_spec  = np.vstack([model_waves, model_flux])
    return model_spec


def plot_tuple_error(intuple, label):
    plt.errorbar(intuple[0], intuple[1], intuple[2], label = label)
    return

def plot_tuple(intuple, label, color, linestyle= 'default'):
    plt.plot(intuple[0], intuple[1], label=label, color= color, linestyle= linestyle)
    return

def get_photon_energy(wavelength):
    return (const.c*const.h/ wavelength).cgs

def convolve_with_passband(input_spec, passband_string):
    passband_dict = {"g":dr2_G_tuple,
                     "rp": dr2_RP_tuple,
                     "bp": dr2_BP_tuple}
    passband_tuple = passband_dict[passband_string]
    input_spec= spt.clean_spectrum(input_spec, np.nanmin(passband_tuple[0]), np.nanmax(passband_tuple[0]), []) #trimming spectrum to wavelengths of the passband
    interpolator = scinterp.CubicSpline(passband_tuple[0], passband_tuple[1])
    interpolated_transmission = interpolator(input_spec[0])
    #plt.plot(input_spec[0], interpolated_transmission, label = passband_string)
    #plt.plot(input_spec[0], input_spec[1]/np.nanmax(input_spec[1]), label = 'data')
    #plt.legend()
    #plt.show()
    #plt.plot(input_spec[0], input_spec[1], label = "Model")
    input_flux = np.copy(input_spec[1])*u.erg/u.s/u.cm**2/u.cm
    input_waves = np.copy(input_spec[0])*u.angstrom
    input_flux= interpolated_transmission*input_flux
    #input_flux= input_flux
    #print ("before anything", input_flux.unit)
    #plt.plot(input_waves, input_flux, label = "transmitted")
    #plt.legend()
    #plt.show()
    delta_lambda = input_waves-np.roll(input_waves,1)
    #input_spec = input_spec[:, 1:]
    input_flux = input_flux[1:]
    input_waves = input_waves[1:]
    delta_lambda= delta_lambda[1:]
    input_flux = input_flux*delta_lambda
    #print("times wavelength", input_flux.unit)
    photon_energies = get_photon_energy(input_waves)
    input_photons = input_flux/photon_energies
    #print("divided by energy:",  input_photons.unit)
    input_photons = input_photons*area #now it's in units of photons/s
    #input_photons = input_photons
    #print("times area:",  input_photons.unit)
    #print delta_lambda
    summed_flux = np.sum(input_photons) #photons/s in the telescope
    #print(passband_string, summed_flux.cgs)
    #return summed_flux
    return summed_flux.si
    
    

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


#model_G_w, model_G_f= convolve_with_passband(model_spec, 'G')
#model_BP_w, model_BP_f= convolve_with_passband(model_spec, 'BP')
#model_RP_w, model_RP_f= convolve_with_passband(model_spec, "RP")

def inverse_distance_modulus(apparent_mag, abs_mag):
    """
    Basically returns the distance required for the given apparent mag for the provided absolute
    magnitude. For use with making radii.
    """
    return (10**((apparent_mag- abs_mag+5)/5.))*u.pc

def get_radius(absolute_mag, teff=teff, logg= logg, passband_string= 'g'):
    #wd=wdatmos.wdmodel(filename='ELM.hdf5')
    #model = wd(Teff = teff, logg = logg)
    ##print("model['flux']:", model['flux'])
    #model_waves = model['w'].data
    #model_flux = model['flux'].data #since we'll be arbitrarily-ish scaling this it won't work.
    #model_spec  = np.vstack
    model_spec = get_model_spec(teff = teff, logg=logg)
    model_band_flux = convolve_with_passband(model_spec,passband_string)
    band_mag= get_mag(model_band_flux.value, passband_string)
    radius = inverse_distance_modulus(band_mag, absolute_mag).to(u.Rsun)
    return radius
    

def get_model_bp_rp(logg=logg, teff=teff, radius = radius):
    #wd=wdatmos.wdmodel(filename='ELM.hdf5')
    #model = wd(Teff = teff, logg = logg)
    #model_waves = model['w'].data
    #model_flux = model['flux'].data #since we'll be arbitrarily-ish scaling this it won't work.
    #model_spec  = np.vstack([model_waves, model_flux])
    model_spec = get_model_spec(logg= logg, teff= teff)
    model_BP =convolve_with_passband(model_spec, 'bp')
    model_RP = convolve_with_passband(model_spec, 'rp')
    BPmag = get_mag(model_BP.value, 'bp')
    RPmag= get_mag(model_RP.value, 'rp')
    BP_RP = BPmag- RPmag
    return BP_RP

def get_model_absmag(logg= logg, teff= teff, radius = radius, passband_string = 'g'):
    #wd=wdatmos.wdmodel(filename='ELM.hdf5')
    #model = wd(Teff = teff, logg = logg)
    #model_waves = model['w'].data
    #model_flux = model['flux'].data 
    ##model_flux = model['flux'].data/4. #apparently it's 4 times the eddington flux
    #model_spec  = np.vstack([model_waves, model_flux])
    model_spec = get_model_spec(logg= logg, teff= teff)
    wavelengths = np.arange(np.nanmin(model_spec[0]), np.nanmax(model_spec[0]), 0.1)
    #fluxes = scinterp.interp1d(wavelengths)
    #fluxes = np.interp(wavelengths, model_spec[0], model_spec[1])
    interpolator= scinterp.CubicSpline(model_spec[0], model_spec[1])
    fluxes = interpolator(wavelengths)
    model_spec = np.vstack([wavelengths, fluxes])
    model_band_flux =convolve_with_passband(model_spec, passband_string)
    mag = get_mag(model_band_flux.value, passband_string)
    absmag= distance_modulus(mag, (radius.to(u.pc)).value)
    return absmag


#def get_BB_luminosity(teff=teff, radius=radius):
    

def get_model_CMD_loc(logg= logg, teff = teff, radius = radius):
    bp_rp = get_model_bp_rp(logg= logg, teff= teff, radius = radius)
    g_absmag= get_model_absmag(logg=logg, teff = teff, radius = radius)
    return g_absmag, bp_rp

def get_mass(radius, logg):
    return 10**(logg)*(u.cm/u.s**2)*radius**2/const.G

def get_radius_from_mass(mass, logg):
    return (np.sqrt(mass*const.G/(10**(logg)* u.cm/u.s**2))).to(u.Rsun)

def calc_phot_bp_rp_excess(g_flux, bp_flux, rp_flux):
    """
    Manual calculation of the 'phot_bp_rp_excess_factor' from the original table
    
    """
    return g_flux/(bp_flux+rp_flux)



######
#plt.axhline(y=0, color ='k')

#plot_tuple(dr2_G_tuple, "G DR2", color = 'g', linestyle = '--')
##plot_tuple(rev_G_tuple, "G _DR2 revised", color = 'g')

#plot_tuple(dr2_RP_tuple, "$G_{RP}$ DR2", color = 'r', linestyle = '--')
##plot_tuple(rev_RP_tuple, "$G _{RP}$DR2 revised", color = 'r')

##plot_tuple(dr2_BP_tuple, "$G_{BP}$ DR2", color = 'b', linestyle = '--')
#plot_tuple(rev_BP_tuple, "$G _{BP}$ DR2 revised", color = 'b')
#plt.plot(model_spec[0], model_spec[1], color = 'k', label = "Teff= " + str(teff) + ", logg= " + str(logg))
#plt.plot(model_G[0], model_G[1], color = 'g', label = "G Teff= " + str(teff) + ", logg= " + str(logg), linestyle = 'none', marker = '.')
#plt.plot(model_BP[0], model_BP[1], color = 'b', label = " BP Teff= " + str(teff) + ", logg= " + str(logg))
#plt.plot(model_RP[0], model_RP[1], color = 'r', label = "RP Teff= " + str(teff) + ", logg= " + str(logg))
#plt.plot(model_spec[0], model_spec[1], color = 'k', label = "Teff= " + str(teff) + ", logg= " + str(logg))
#plt.plot(model_G_w, model_G_f, color = 'g', label = "G Teff= " + str(teff) + ", logg= " + str(logg))
#plt.plot(model_BP_w, model_BP_f, color = 'b', label = " BP Teff= " + str(teff) + ", logg= " + str(logg))
#plt.plot(model_RP_w, model_RP_f, color = 'r', label = "RP Teff= " + str(teff) + ", logg= " + str(logg))
#plt.xlabel(r'Wavelength ($\AA$)')
#plt.ylabel('Flux')
#plt.legend()
#plt.show()
