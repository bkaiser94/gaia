"""
Created by Ben Kaiser (UNC-Chapel Hill) 2023-06-05.

This should handle the Gaia XP spectra in their weird Hermite polynomial coefficient form and
make it something somewhat palatable. Whether that is a sampled version or just at least some
means of plotting it in pyplot. Just something so I can at least look at the data.

It also seems decently likely that this will handle the retrieval of the spectra too since they
aren't included with the general retrieval by default and have to be retrieved via datalink.





"""

from __future__ import print_function


#import matplotlib
#matplotlib.use('pdf')


import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, QTable
import matplotlib.pyplot as plt
import scipy.stats as scistats
#import seaborn as sns
#import astropy
import gaiaxpy as xpy
import pandas as pd
import sys


#import passband_model_convolution as pmc
#import gaia_extinction
#import wdatmos
import plotting_dicts as pod

output_dir='GaiaDR3_XP_spectra/WDMS_wide_binaries/'
#input_file='20210201_DZs_for_J1636paper_gaia_gaiaDR3.csv'
#input_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf.csv'
input_file='WDJ1948m1011_gaiaDR3.csv'
#input_file='LTT3218_GaiaDR3.csv'
credentials_file= 'Gaia_credentials.txt'

#output_file='LHS2534_GaiaDR3_XPspectrum.csv'
output_file='wdbinary_GaiaDR3_XPspectrum.csv'
save_file=False
save_plots=False
single_index=3

credentials_frame=np.genfromtxt(credentials_file,dtype=str)
input_table=Table.read(input_file)
input_table.pprint()
has_xp=np.where((np.array(input_table['has_xp_continuous'])=='True' )|( np.array(input_table['has_xp_continuous'])=='true'))
sub_table=input_table[has_xp]
print('sub_table')
sub_table.pprint()

source_list=[]
name_list=[]
for row in sub_table:
    try:
        source=row['DESIGNATION'].split(' ')[-1]
    except KeyError:
        source=row['designation'].split(' ')[-1]
    source_list.append(source)
    name_list.append(row['name'])
    

#if len(sub_table)>1:
    #for row in sub_table:
        ##source=row['DESIGNATION'].split(' ')[-1]
        #source=row['designation'].split(' ')[-1]
        #source_list.append(source)
        #name_list.append(row['name'])
#else:
    #source=sub_table['designation'][0].split(' ')[-1]
    #source_list.append(source)
    #name_list.append(sub_table['name'][0])

print(source_list)
print(credentials_frame.shape)
print('about to retrieve and calibrate spectra')
calibrated_spectra, sampling=xpy.calibrate(source_list,username=credentials_frame[0],password=credentials_frame[1],truncation=True)
print('spectra retrieved and calibrated')

print(calibrated_spectra)

print('sampling',sampling)

conversion_factor=1e18 #multiple by which you have to multiply the flux in W/m^2/nm/s to get to erg/cm^2/Anstrom/s *10^-16 as used in my Goodman spectra
#for row, name in zip(calibrated_spectra,name_list):
    #row['source_id']=name
calibrated_spectra['source_id']=name_list
xpy.plot_spectra(calibrated_spectra,sampling=sampling,multi=True) #this plots all of the spectra at one time.

#print(calibrated_spectra['flux'].units)
#plt.plot(sampling,calibrated_spectra)
calibrated_spectra['flux']=calibrated_spectra['flux']*conversion_factor
calibrated_spectra['flux_error']=calibrated_spectra['flux_error']*conversion_factor
sampling=sampling*10.
xpy.plot_spectra(calibrated_spectra,sampling=sampling,multi=True) #this plots all of the spectra at one time.
print(sampling.shape)
print(calibrated_spectra['flux'].shape)
plt.xlabel('Wavelength (Angstroms')
plt.ylabel('Flux (cgs units)')
plt.show()

#sys.exit()
print('type(calibrated_spectra)',type(calibrated_spectra))
for index,calibrated_spectrum in calibrated_spectra.iterrows():
    print('calibrated_spectrum', calibrated_spectrum)
    ra_dec=coord.SkyCoord(ra=sub_table[index]['ra'], dec=sub_table[index]['dec'], unit=(u.deg,u.deg),frame='icrs')
    ra_dec_string=ra_dec.to_string(style='hmsdms')
    full_radec=ra_dec.to_string(style='hmsdms')
    replace_chars=['h','m','d']
    for char in replace_chars:
        full_radec=full_radec.replace(char,':')
    full_radec=full_radec.replace('s','')
    ra_dec_list=full_radec.split(' ')
    ra_string=ra_dec_list[0][:5]
    dec_string=ra_dec_list[1][:6]
    ra_string=ra_string.replace(':','')
    dec_string=dec_string.replace(':','')
    #ra_string.replace('s','')
    #dec_string.replace('s','')
    #dec_string=dec_string.replace('-','m')
    #dec_string=dec_string.replace('+','p')

    #plt.title(name_list[index]+' '+calibrated_spectra['source_id']+', '+ra_dec.to_string(style='hmsdms')+', ' +'G='+str(sub_table[index]['phot_g_mean_mag'])[:5])
    if save_plots:
        plt.figure(figsize=(1036./540*6.,6.))
    else:
        pass
    #plt.title('index:'+str(index)+', '+name_list[index]+' '+str(sub_table[index]['source_id'])+', '+ra_dec.to_string(style='hmsdms')+', ' +'G='+str(sub_table[index]['phot_g_mean_mag'])[:5])
    #plt.title('index:'+str(index)+', '+'J'+ra_string+dec_string+', '+str(sub_table[index]['source_id'])+', '+full_radec+', ' +'G='+str(sub_table[index]['phot_g_mean_mag'])[:5])
    plt.title(calibrated_spectrum['source_id'])
    plt.plot(sampling, calibrated_spectra['flux'][index])
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Flux (cgs units)')
    #if sub_table[index]['phot_g_mean_mag']<18.:
        #plt.plot(sampling, calibrated_spectra['flux'][index])
        #plt.xlabel('Wavelength (Angstroms')
        #plt.ylabel('Flux (cgs units)')
        #plt.errorbar(sampling, calibrated_spectra['flux'][index],yerr=calibrated_spectra['flux_error'][index])
    #xpy.plot_spectra(calibrated_spectrum,sampling=sampling,multi=False)
        #plt.show() 
    
    if save_plots:
        print('saving plot')
        for char in replace_chars:
            ra_dec_string=ra_dec_string.replace(char,'')
        ra_dec_list=ra_dec_string.split(' ')
        ra_string=ra_dec_list[0][:4]
        dec_string=ra_dec_list[1][:5]
        ra_string.replace('s','')
        dec_string.replace('s','')
        dec_string=dec_string.replace('-','m')
        dec_string=dec_string.replace('+','p')
        plt.savefig(output_dir+'J'+ra_string+dec_string+'_spec'+str(index)+'_wide_binary_DR3XP.pdf')
        plt.clf()
    else:
        print('showing plot')
        plt.show()


print(calibrated_spectra['flux'][0])

single_spec_flux=calibrated_spectra['flux'][single_index]
single_spec_error=calibrated_spectra['flux_error'][single_index]*conversion_factor
#single_spec_flux=single_spec_flux.to(u.erg/(u.cm**2)/u.angstrom/u.s)
#single_spec_flux=single_spec_flux.cgs
print(single_spec_flux)
print(calibrated_spectra['source_id'][single_index])

wavelengths=sampling*10. #default sampling from Gaia is in nm and I want Angstroms.

output_spec=np.vstack([wavelengths,single_spec_flux,single_spec_error])

if save_file:
    print('saving output spectrum',output_file)
    np.savetxt(output_dir+output_file,output_spec)
    print('file saved.')
else:
    pass
    








