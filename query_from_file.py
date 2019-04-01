"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-12-31

This should specifically be capable of taking in  a list of coordinates and such and then querying the gaia archive for the associated objects... hopefully. We'll see! This has only failed every other time I've attempted it.

"""



import numpy as np
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column
 
 
filter_confused_sources= True
search_radius = 4.08 #in arcseconds
#search_radius =  10.
 
 
#input_file= 'BLAPs_table1_Piet2017.csv'
#output_file= 'BLAPs_gaia.csv'
#input_file= 'pulsar_companions.txt'
#output_file= 'pulsar_companions_gaia.csv'
#input_file= 'hot_wind_wds.txt'
#output_file= 'hot_wind_wds_gaia.csv'

#input_file= 'hv_wds.txt'
#output_file= 'hv_wds_gaia.csv'

#input_file= 'pre_elms.txt'
#output_file= 'pre_elms_gaia.csv'

#input_file= 'elm_survey.txt'
#input_file='apj522588t5_mrt.txt'
#output_file= 'elm_survey_gaia.csv'

#input_file= '20190104_chris.csv'
input_file= 'exc1_8_2_2_purple_search_gmaglimit.csv'
input_file= 'usdMs.csv'
input_file='Lindegren_odd_survivors.csv'
#input_file= 'dC_sample_roulston2018.csv'
#input_file='20190111_red_things.csv'
#input_file= 'alt_red_things_g_rp_greater_17.csv'
#input_file='WD_cooling_tip.csv'
#input_file= 'weird_CPM_binary.csv'
#input_file='Lindegren_appC_selB_antiC_nobulgedisk.csv'
#input_file='Lindegren_appC_selB_antiC_cut2.csv'
#input_file='20190109_blue.csv'
#input_file='20190121_excess_interesting.csv'
#input_file='ar_sco.txt'
#input_file='20190123_new_red_things.csv'
#input_file= 'DQpec.csv'
#input_file= '20190128_blue_line.csv'
#input_file='20190128_wdMS.csv'
#input_file= 'SLBs.csv'
#input_file = 'RNe.csv'
#input_file= "FA_ultracool_WD_attempts.csv"
#input_file= 'coolDZ_Na.csv'
#input_file='two_low_reds.csv'
#input_file= 'mansergas_2010_pceb.csv'
#input_file= 'esdM_cand_subset.csv'
#input_file= 'Gianinas2016_coolWDs.txt'
#input_file='Eriks_disk_candidates.csv'
#input_file='all_l-0.3bp_g_gaia_corr.csv'

if filter_confused_sources:
    output_name_parts = input_file.split('.')
    output_file= output_name_parts[0]+ '_gaia_sc.' + output_name_parts[1]
else:
    output_name_parts = input_file.split('.')
    if output_name_parts[1]=='txt':
        output_name_parts[1]= 'csv'
    else:
        pass
    output_file= output_name_parts[0]+ '_gaia.' + output_name_parts[1]
#output_file='l-0.3bp_g_gaia_corr_full.csv'
credentials_file= 'Gaia_credentials.txt'
print("Logging in.")
Gaia.login(credentials_file= credentials_file)
print("Log in succesful.")


coord_list= []
####allarray= np.genfromtxt(input_file, names=True, delimiter=',', dtype='U6')
allarray= np.genfromtxt(input_file, delimiter= ',', dtype= str)
#allarray=np.genfromtxt(input_file, dtype =str, delimiter= ' \t',skip_header=2)
#####allarray= np.genfromtxt(input_file, delimiter=[11,13,13], dtype= str)
####allarray= np.loadtxt(input_file, delimiter=',', dtype=bytes).astype(str)

#print(allarray)
#ra_array = allarray['Ra']
#dec_array= allarray['Dec'] #Now I've got the RAs and Decs
allarray=allarray[1:]
#print(allarray)

ra_array=allarray[:,1]
dec_array= allarray[:,2] #I couldn't get the names and dtype to work so I'm just doing it the less efficient way I guess
name_array=allarray[:,0]

collected_results=[]
print(name_array)
print(ra_array)
print(dec_array)

def cone_search(ra, dec):
    #print(ra)
    #print(dec)
    print(ra)
    if ((":" in ra) or ( " " in ra) or ( "\t" in ra)):
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
    j= Gaia.cone_search_async(coordinate, radius)
    r=j.get_results()
    r.pprint()
    return r

for ra,dec,name in zip(ra_array, dec_array,name_array):
    results= cone_search(ra,dec)
    #print(results.shape)
    #new_col= Column(name, name='name')
    #print(results[0])
    #print(results[0])
    #output_row= results[0]
    name=str(name)
    #name='SDSS'+name
    print(name)
    print(type(name))
    try:
        table_length = len(results['dist'])
        print(table_length)
        if (filter_confused_sources and table_length>1):
            print("Too many sources in aperture")
        else:
            #new_array = np.full(table_length, name, dtype=str)
            #new_array= np.empty(table_length, dtype=str)
            #new_array[:]=name
            #print(new_array)
            new_array=[]
            for i in range(0,table_length):
                new_array.append(name)
            print(new_array)
            name_col= Column(new_array, name='name',dtype=str) #yeah that's a confusing series of 'names'
            #output_row.add_column(name_col)
            results.add_column(name_col)
            try:
                collected_results.append(results[0])
                #collected_results.append(output_row)
                #collected_results.append(results)
            except IndexError:
                print("index error")
                pass
    except TypeError as error:
        print(error)
        pass
    
stacked_results= vstack(collected_results)
stacked_results.pprint()

stacked_results.write(output_file, format='ascii.csv', overwrite=True)

#for thing in coord_list:
    #print(thing.ra.to(u.hourangle), thing.dec)
