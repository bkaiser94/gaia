"""
Created by Ben Kaiser (UNC-Chapel Hill) 2023-06-05.

This should handle the Gaia XP spectra in their weird Hermite polynomial coefficient form and
make it something somewhat palatable. Whether that is a sampled version or just at least some
means of plotting it in pyplot. Just something so I can at least look at the data.

It also seems decently likely that this will handle the retrieval of the spectra too since they
aren't included with the general retrieval by default and have to be retrieved via datalink.





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
import gaiaxpy as xpy
import pandas as pd


#import passband_model_convolution as pmc
#import gaia_extinction
#import wdatmos
import plotting_dicts as pod

output_dir='GaiaDR3_XP_spectra/'
input_file='20210201_DZs_for_J1636paper_gaia_gaiaDR3.csv'
#input_file='LTT3218_GaiaDR3.csv'
credentials_file= 'Gaia_credentials.txt'
output_file='LHS2534_GaiaDR3_XPspectrum.csv'
save_file=False
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
print(credentials_frame[0])
print(credentials_frame[1])
print('about to retrieve and calibrate spectra')
calibrated_spectra, sampling=xpy.calibrate(source_list,username=credentials_frame[0],password=credentials_frame[1],truncation=True)
print('spectra retrieved and calibrated')

print(calibrated_spectra)

print(sampling)

conversion_factor=1e18 #multiple by which you have to multiply the flux in W/m^2/nm/s to get to erg/cm^2/Anstrom/s *10^-16 as used in my Goodman spectra
#for row, name in zip(calibrated_spectra,name_list):
    #row['source_id']=name
calibrated_spectra['source_id']=name_list
#print(calibrated_spectra['flux'].units)
#plt.plot(sampling,calibrated_spectra)
calibrated_spectra['flux']=calibrated_spectra['flux']*conversion_factor
xpy.plot_spectra(calibrated_spectra,sampling=sampling,multi=True)
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
    








