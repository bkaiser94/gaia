"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-12-31

This should specifically be capable of taking in  a list of coordinates and such and then querying the gaia archive for the associated objects... hopefully. We'll see! This has only failed every other time I've attempted it.

"""



import numpy as np
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack
 
 
input_file= 'BLAPs_table1_Piet2017.csv'
output_file= 'BLAPs_gaia.csv'


credentials_file= 'Gaia_credentials.txt'
print("Logging in.")
Gaia.login(credentials_file= credentials_file)
print("Log in succesful.")


#allarray= np.genfromtxt(input_file, names=True, delimiter=',', dtype='U6')
allarray= np.genfromtxt(input_file, delimiter= ',', dtype= str)
#allarray= np.loadtxt(input_file, delimiter=',', dtype=bytes).astype(str)

#print(allarray)
#ra_array = allarray['Ra']
#dec_array= allarray['Dec'] #Now I've got the RAs and Decs
allarray=allarray[1:]
#print(allarray)

ra_array=allarray[:,1]
dec_array= allarray[:,2] #I couldn't get the names and dtype to work so I'm just doing it the less efficient way I guess

collected_results=[]
print(ra_array)

def cone_search(ra, dec):
    #print(ra)
    #print(dec)
    coordinate = coord.SkyCoord(ra = ra, dec =dec, unit = (u.hourangle, u.deg), frame = 'icrs')
    #coordinate = coord.SkyCoord(ra = ra*u.hourangle, dec =dec*u.degree, frame = 'icrs')
    radius = 1.5*u.arcsecond
    j= Gaia.cone_search_async(coordinate, radius)
    r=j.get_results()
    r.pprint()
    return r


for ra,dec in zip(ra_array, dec_array):
    results= cone_search(ra,dec)
    #print(results[0])
    #print(results[0])
    try:
        collected_results.append(results[0])
    except IndexError:
        pass
    
stacked_results= vstack(collected_results)
stacked_results.pprint()

stacked_results.write(output_file, format='ascii.csv', overwrite=True)
