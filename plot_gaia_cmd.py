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

num_targs = 100000
num_targs = int(num_targs)
generic_input = 'top'+str(num_targs) + '_nearby_gaia.csv'
target_input = 'target_gaia.csv'

#generic_table = Table.read('top500_nearby_gaia.csv')
#generic_table = Table.read('top5000_nearby_gaia.csv')
generic_table = Table.read(generic_input)
target_table = Table.read(target_input)



#generic_table.pprint()


generic_parallax = generic_table ['parallax']
generic_parallax = generic_parallax *1e-3 #parallax in arcseconds now
generic_distance = 1./generic_parallax #parsec distance

generic_extinction = generic_table['a_g_val']

generic_g_mag = generic_table['phot_g_mean_mag']
generic_bp_rp = generic_table['bp_rp']

#print(target_table['ra'])
#print(target_table['dec'])
target_parallax = target_table['parallax']
target_parallax = target_parallax*1e-3
target_distance = 1./target_parallax

target_extinction= target_table['a_g_val']

target_g_mag = target_table['phot_g_mean_mag']
target_bp_rp = target_table['bp_rp']



def distance_modulus(g_mag, distance, extinction = 0.0):
    return g_mag - 5*np.log10(distance/10.)
    #return g_mag - 5*np.log10(distance/10.)- np.float_(extinction)


generic_g_absmag = distance_modulus(generic_g_mag, generic_distance, extinction = generic_extinction)
target_g_absmag = distance_modulus(target_g_mag, target_distance,extinction= target_extinction)


#print(generic_parallax)
print(target_g_mag)
print(target_extinction)
print(target_table['a_g_percentile_lower'],target_table['a_g_percentile_upper'])
print(target_g_absmag)
print(target_bp_rp)

#plt.plot(generic_parallax, linestyle = 'none', marker = 'o')
#plt.show()

#plt.hist(generic_parallax)
#plt.show()

#plt.hist(generic_distance)
#plt.xlabel('Distance (pc)')
#plt.show()

def make_density_plot(g_abs, bp_rp):
    #Calculate the point density
    xy = np.vstack([np.array(bp_rp),np.array(g_abs)])
    z = scistats.gaussian_kde(xy)(xy)

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
plt.plot(target_bp_rp, target_g_absmag, marker = '*', markersize = 8, color = 'b')
#plt.hist2d(generic_bp_rp, generic_g_absmag, bins = 100, cmap = 'Reds')
plt.ylim([-4, 16])

plt.gca().invert_yaxis()
plt.xlabel(r'$G_{BP} - G_{RP}$')
plt.xlim([-1,5])
plt.ylabel(r'$M_G$')
plt.show()

#plt.hist2d(generic_bp_rp, generic_g_absmag, bins = 100, cmap = 'Reds')
#plt.gca().invert_yaxis()
#plt.xlabel(r'$G_{BP} - G_{RP}$')
#plt.ylabel(r'$G$')
#plt.show()
