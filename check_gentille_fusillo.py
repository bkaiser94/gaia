"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-05-17

So I changed my mind. This actually just does a cross-match with Nicola's Gaiawd catalogue directly in Gaia
coordinates to see if there's a common object, and then retains a field that indicates it was cross-matched, the Pwd value, and aspirationally the f_Pwd flag, but for whatever reason the vizier astroquery tool can't handle that at all.

"""

import numpy as np
#from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column
import time
start = time.time()

#dictionary of acceptable values and associated acceptability cuts
#astrometric_sigma5d_max < 1.5

input_file= '20190516B_retargeted_purple_search_gaia_scbd.csv'
input_table = Table.read(input_file)

search_radius = 3#in arcseconds


cat_name = 'J/MNRAS/482/4570/gaia2wd'

ra_array = input_table['ra']
dec_array = input_table['dec']
name_array = input_table['name']
num_targets= ra_array.shape[0]

new_cols = ['gf19',
            'pwd',
            'f_pwd',
            'sdss_spec']

for new_col in new_cols:
    try:
        print('Test repeat column: ', input_table[new_col][0])
    except KeyError:
        print(new_col, " column doesn't exist, so we're making it.")
        new_length = len(input_table['dist'])
        if new_col=='gf19':
            false_array = np.full((new_length,), False, dtype=bool)
            new_column = Column(false_array, name=new_col,dtype=bool)
        elif new_col=='pwd':
            new_array= np.full((new_length,), 0, dtype=float)
            new_column=Column(new_array, name=new_col,dtype=float)
        elif new_col == 'f_pwd':
            new_array= np.full((new_length,), 0, dtype=int)
            new_column=Column(new_array, name=new_col,dtype=int)
        elif new_col == 'sdss_spec':
            new_array= np.full((new_length,), '', dtype=str)
            new_column=Column(new_array, name=new_col,dtype=str)
        else:
            print("no matching names:", new_col)
        input_table.add_column(new_column)

def find_time_difference(begin, end):
    difference = end-begin
    hours = int(difference) // 3600
    minutes = int(difference%3600)//60
    seconds = difference%60
    #print("Runtime: ", hours, 'h', minutes, 'm', seconds, 's')
    output_string= ' '.join(["Runtime: ", str(hours), 'h',str( minutes), 'm', str(seconds), 's'])
    return output_string


def cone_search(ra, dec, search_radius= search_radius):
    #print(ra)
    #print(dec)
    print(ra)
    if ((":" in str(ra)) or ( " " in str(ra)) or ( "\t" in str(ra))):
        print("RA in hour angle")
        coordinate = coord.SkyCoord(ra = ra, dec =dec, unit = (u.hourangle, u.deg), frame = 'icrs')
    else:
        coordinate = coord.SkyCoord(ra = ra, dec =dec, unit = (u.deg, u.deg), frame = 'icrs')
    #coordinate = coord.SkyCoord(ra = ra*u.hourangle, dec =dec*u.degree, frame = 'icrs')
    #radius = 1.5*u.arcsecond
    radius= search_radius*u.arcsecond
    #coord_list.append(coordinate)
    #print(coordinate.ra.to(u.hourangle), coordinate.dec)
    #radius = 5*u.arcsecond
    #j= Gaia.cone_search_async(coordinate, radius)
    j= Vizier.query_region(coordinate, radius=radius, catalog=cat_name)
    #r=j.get_results()
    #j.pprint()
    print(j)
    matched_result= False
    try:
        print(j[0])
        for thing in j[0]:
            print(thing)
        matched_result=True
        print(j[0]['Pwd'])
        pwd= j[0]['Pwd']
        sp= j[0]['Sp']
        try:
            f_pwd= j[0]['f_Pwd']
        except KeyError as error:
            print("KeyError:", error)
            f_pwd= 0
    except IndexError:
        pwd=0
        f_pwd=0
        sp=''
        pass
    return matched_result,pwd, f_pwd,sp

#GF19_matched= 
target_num = 1
for ra,dec,name, ogf19, opwd, of_pwd in zip(ra_array, dec_array,name_array, input_table['gf19'], input_table['pwd'], input_table['f_pwd']):
    begin_time = time.time()
    print("\n\ntarget ", str(target_num) + '/' + str(num_targets),'\n\n')
    #results= cone_search(ra,dec)
    gf19, pwd, f_pwd, sp= cone_search(ra,dec)
    #print(results.shape)
    #new_col= Column(name, name='name')
    #print(results[0])
    #print(results[0])
    #output_row= results[0]
    name=str(name)
    #name='SDSS'+name
    print(name)
    print(type(name))
    print('pwd', pwd)
    row= input_table[target_num-1]
    row['pwd']=pwd
    row['f_pwd']=f_pwd
    row['gf19']=gf19
    print('sp', sp)
    print('type',type(sp))
    #print('dtype',sp.dtype)
    try:
        row['sdss_spec']= sp
    except ValueError as error:
        print("ValueError:", error)
        row['sdss_spec']=''
    #try:
        ##table_length = len(results['dist'])
        #print('trying table length')
        #table_length= len(results['ra'])
        #print(table_length)
        #if (filter_confused_sources and table_length>1):
            #print("Too many sources in aperture")
        #else:
            #if (filter_bright_contam or filter_dense):
                #bright_results= cone_search(ra,dec, search_radius=bright_search_radius)
                #if (filter_bright_contam and np.nanmin(bright_results['phot_g_mean_mag']) < bright_star_limit):
                    #print("Star brighter than G=" + str(bright_star_limit) + " , so likely contamination, thus excluding.")
                #elif (filter_dense and len(bright_results['dist']) >= dense_number):
                    #print("Stellar density greater than 100,000 deg^-2, there are ", len(bright_results['dist']), " stars")
                #else:
                    ##new_array = np.full(table_length, name, dtype=str)
                    ##new_array= np.empty(table_length, dtype=str)
                    ##new_array[:]=name
                    ##print(new_array)
                    #new_array=[]
                    #for i in range(0,table_length):
                        #new_array.append(name)
                    #print(new_array)
                    #name_col= Column(new_array, name='name',dtype=str) #yeah that's a confusing series of 'names'
                    ##output_row.add_column(name_col)
                    #results.add_column(name_col)
                    #try:
                        #collected_results.append(results[0])
                        ##collected_results.append(output_row)
                        ##collected_results.append(results)
                    #except IndexError:
                        #print("index error")
                        #pass
            #else:
                ##new_array = np.full(table_length, name, dtype=str)
                ##new_array= np.empty(table_length, dtype=str)
                ##new_array[:]=name
                ##print(new_array)
                #new_array=[]
                #for i in range(0,table_length):
                    #new_array.append(name)
                #print(new_array)
                #name_col= Column(new_array, name='name',dtype=str) #yeah that's a confusing series of 'names'
                ##output_row.add_column(name_col)
                #results.add_column(name_col)
                #try:
                    #collected_results.append(results[0])
                    ##collected_results.append(output_row)
                    ##collected_results.append(results)
                #except IndexError:
                    #print("index error")
                    #pass
    #except TypeError as error:
        #print('TypeError:', error)
        #pass
    target_num += 1
    end_time= time.time()
    print('Target search time: ', find_time_difference(begin_time, end_time))
    print('\nTotal elapsed time', find_time_difference(start, end_time))
    print('\n\n\n===============')
    
    
#stacked_results= vstack(collected_results)
#stacked_results.pprint()
input_table.pprint()
input_table['gf19'].pprint()

input_table.write(input_file, format='ascii.csv', overwrite=True)
#stacked_results.write(output_file, format='ascii.csv', overwrite=True)

end_overall= time.time()
print('\n\n')
print('Total elapsed time: ', find_time_difference(start, end_overall))
