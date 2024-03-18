"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-11-08 from query_from_file.py

This should take a Gaia DR2 file and then query DR3 for the updated information, mainly with 
the goal of using this for kinematics, but I'm going to have it grab everything because why not. 
It will be using RA and DEC from DR2 to query in DR3, which means it will be using 2015.5 
epoch coordinates on 2016 epoch coordinates. I don't have any immediate intention of moving 
them forward in time to match because I figure I might as well get the Gaia ID cross matching 
working at that point. I should probably be doing that anyway, but time-crunch so we're going 
anyway.

"""



import numpy as np
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column
import matplotlib.pyplot as plt

import time
import sys
start = time.time()

#Setting which gaia table is being used for queries.

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" 
#Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"
#Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source"

#tables=Gaia.load_tables(only_names=True)
#for table in tables:
    #print(table.get_qualified_name())


filter_confused_sources= False
#search_radius=1.
search_radius = 4.08 #in arcseconds
#search_radius =  10.

filter_bright_contam= False
filter_dense= False#filter to remove crowded fields; didn't use 'crowded' because it has a 'c' like confusion
 
filter_sepbright_func=False
sepbright_radius=70. #radius in arcseconds out to which the query should retrieve objects. This is going to be slow.
 
#filter_confused_sources= False
#search_radius=1.
##search_radius = 4.08 #in arcseconds
##search_radius =  10.

#filter_bright_contam= False
#filter_dense= False#filter to remove crowded fields; didn't use 'crowded' because it has a 'c' like confusion

max_stellar_density=110000. #number of stars per square degree above which a field is considered to crowded if filter_dense=True
 
bright_search_radius = 15.
bright_star_limit= 15 #faintest G mag of what's considered a bright star for purposes of filtering
#dense_number= 6 #corresponds to >110,008 stars per square degree for 15" search radius, so should be greater than or equal to this number to qualify as 'dense' or 'crowded'#This is now calculated in the code using find_star_limit(bright_search_radius), using the max_stellar_density as the filtering method.



#################### binary separation function values#########

close_mag_diff=-4.5 #constant value inside the closest radius (these are basically completely incapable of being disentangled in the reduction (or rather removing the contaminant flux so the contaminant flux needs to be much less than the target. Hence the negative value
close_mag_radius=4. #in arcseconds
mid_mag_diff=4.6 
mid_mag_out_radius=13. #in arcseconds

mag_coeff=-9.8
mag_center=7.
mag_spread=15.
mag_offset=14.6




#################################################


#name_index=0
#ra_index=1
#dec_index=2

name_index=2
ra_index=0
dec_index=1



 
 
 
#input_file= 'gas_disk_objs.csv'
#input_file='Kilic_velocities_tabsupp_2_rand50.csv'
#input_file='20210305B_ultracool_switchback_gaia.csv'
#input_file='20220929_reallydim_WD_search.csv'
#input_file='20210201_DZs_for_J1636paper_gaia.csv'
#input_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWDradec.csv'
#input_file='20210201_DZs_for_J1636paper_gaia_gaiaDR3.csv'

#input_file='20191024_DZNa_and_redder_search.csv'
#input_file='Torres_50_cat_original.csv'
#input_file= '20190516B_retargeted_purple_search.csv'
#input_file= 'WDpec.csv'
#input_file= '20190829_alkaliWD_targeted.csv'
#input_file='20190917_alkaliWD_attempt2.csv'
#input_file= 'pre_elms.txt'
#output_file= 'pre_elms_gaia.csv'
#input_file='NLTT5306_comps.csv'
#input_file= 'elm_survey.txt'
#input_file='apj522588t5_mrt.txt'
#output_file= 'elm_survey_gaia.csv'
#input_file='20190511_observed_objects.txt'
#input_file= '20190104_chris.csv'
#input_file='20190422_obs_objects.txt'
#input_file='sdssj1330p6435.csv'
#input_file ='sdssj1330_similar.csv'
#input_file ='sdssj1330_similar_subset.csv'
#input_file= 'exc1_8_2_2_purple_search_gmaglimit.csv'
#input_file= 'usdMs.csv'
#input_file='Lindegren_odd_survivors.csv'
#input_file= 'mdwarf_spTs.csv'
#input_file='20190405_purple_search.csv'
#input_file= 'Lindegren_appC_selC_noBLMCbperr_wd1401_bincomp.csv'
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
#input_file= '20190516_targeted_purple_search.csv'
#input_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_loosecutsonly_justWDradec.csv'
input_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_MSmatches.csv'

additional_suffixes= ''

if filter_confused_sources:
    additional_suffixes='sc'
if filter_bright_contam:
    additional_suffixes= additional_suffixes+'b'
if filter_dense:
    additional_suffixes= additional_suffixes+'d'
if filter_sepbright_func:
    additional_suffixes=additional_suffixes+'_sbf'
else:
    pass
    
if Gaia.MAIN_GAIA_TABLE=='gaiadr3.gaia_source':
    gaia_string='gaiaDR3'
elif Gaia.MAIN_GAIA_TABLE=='gaiaedr3.gaia_source':
    gaia_string='gaiaeDR3'
elif Gaia.MAIN_GAIA_TABLE=='gaiadr2.gaia_source':
    gaia_string='gaiaDR2'
else:
    pass

#if filter_confused_sources:
    #output_name_parts = input_file.split('.')
    ##output_file= output_name_parts[0]+ '_gaia_sc.' + output_name_parts[1]
    ##output_file= output_name_parts[0]+ '_gaiaDR3_sc' + additional_suffixes+'.'+output_name_parts
    #output_file= output_name_parts[0]+ '_'+ gaia_string+ '_sc' + additional_suffixes+'.'+output_name_parts[1]
#else:
    #output_name_parts = input_file.split('.')
    #if output_name_parts[1]=='txt':
        #output_name_parts[1]= 'csv'
    #else:
        #pass
    ##output_file= output_name_parts[0]+ '_gaiaDR3.' + output_name_parts[1]
    #output_file= output_name_parts[0]+ '_'+gaia_string+'.' + output_name_parts[1]

output_name_parts = input_file.split('.')
if output_name_parts[1]=='txt':
    output_name_parts[1]= 'csv'
output_file= output_name_parts[0]+ '_'+ gaia_string+'_'+ additional_suffixes+'.'+output_name_parts[1]
print('\n\noutput_file:',output_file,'\n\n')

#output_file='l-0.3bp_g_gaia_corr_full.csv'
credentials_file= 'Gaia_credentials.txt'
print("Logging in.")
Gaia.login(credentials_file= credentials_file)
print("Log in succesful.")


coord_list= []
####allarray= np.genfromtxt(input_file, names=True, delimiter=',', dtype='U6')
#allarray= np.genfromtxt(input_file, delimiter= ',', dtype= str)

##allarray=np.genfromtxt(input_file, dtype=np.str_, delimiter= ' \t', skip_header=1)
##allarray=np.genfromtxt(input_file, dtype=np.str_, delimiter= '\t')
######allarray= np.genfromtxt(input_file, delimiter=[11,13,13], dtype= str)
#####allarray= np.loadtxt(input_file, delimiter=',', dtype=bytes).astype(str)
#print(allarray)
#print(allarray.shape)
##print(allarray)
##ra_array = allarray['Ra']
##dec_array= allarray['Dec'] #Now I've got the RAs and Decs


#allarray=allarray[1:] #remove the header of the file (or column names)

#print(allarray)

#ra_array=allarray[:,1]
#dec_array= allarray[:,2] #I couldn't get the names and dtype to work so I'm just doing it the less efficient way I guess
#name_array=allarray[:,0]



#ra_array=allarray[:,ra_index]
#dec_array= allarray[:,dec_index] #I couldn't get the names and dtype to work so I'm just doing it the less efficient way I guess
#name_array=allarray[:,name_index]

######### The stuff above was from the original query file, but I'm going to read-in astropy tables, so I don't want numpy arrays.

all_table=Table.read(input_file, delimiter=',')
ra_array=all_table['ra']
dec_array=all_table['dec']
#name_array=all_table['name']


collected_results=[]
print(name_array)
print(ra_array)
print(dec_array)

num_targets= ra_array.shape[0]
print('\n\n')
print('num_targets:', num_targets)
print('\n\n')

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
    j= Gaia.cone_search_async(coordinate, radius)
    r=j.get_results()
    r.pprint()
    return r


def find_star_limit(radius):
    """
    radius: search radius in arc seconds.
    
    will use max_stellar_density value to obtain the maximum number of stars that can be inside the search radius too.
    returned value is the first number of stars for which the field is too crowded. i.e. >= star_limit means too crowded.
    
    """
    radius_deg = (radius*u.arcsecond).to(u.deg).value
    radius_rad=radius_deg*np.pi/180.
    solid_angle=4*np.pi*np.sin(radius_rad/2.)**2. #solid angle in steradians
    solid_angle_sqdeg=solid_angle*((180./np.pi)**2.) #solid angle in square degrees
    star_limit=max_stellar_density*solid_angle_sqdeg
    return star_limit

def far_mag_func(dist_column):
    #print('dist_column',dist_column)
    mag_boundary=mag_coeff*np.exp(-1.*(dist_column-mag_center)**2./(2*mag_spread**2.))+mag_offset
    #print('mag_boundary calculated.')
    return mag_boundary


def check_distance_mag(results):
    """
    results: the results from a cone search, so this contains all of the stars in a field.
    
    
    Needs to check the relative magnitude and distance for each of the 3 layers of the 
    """
    good_field=True #this boolean is what will be returned by the function and should be used as an escape method at the end.
    target_row=results[0]
    other_results=results[1:]
    mag_diff=target_row['phot_g_mean_mag']-other_results['phot_g_mean_mag']
    close_pairs=np.where(other_results['dist']<close_mag_radius)
    close_bad_inds=np.where(mag_diff[close_pairs]>close_mag_diff)
    #plt.scatter(other_results[close_pairs]['dist'],np.ones(other_results[close_pairs]['dist'].shape)*mid_mag_diff,label='close boundary')
    print('close_bad_inds',close_bad_inds)
    if close_bad_inds[0].shape[0]>0:
        print('bad close stars')
        good_field=False
    if good_field:
        mid_pairs=np.where((other_results['dist']>=close_mag_radius) & (other_results['dist'] <= mid_mag_out_radius))
        mid_bad_inds=np.where(mag_diff[mid_pairs]>mid_mag_diff)
        #plt.scatter(other_results[mid_pairs]['dist'],np.ones(other_results[mid_pairs]['dist'].shape)*mid_mag_diff,label='mid boundary')
        print('mid_bad_inds',mid_bad_inds)
        if mid_bad_inds[0].shape[0]>0:
            print('bad mid stars')
            good_field=False
    if good_field:
        far_pairs=np.where(other_results['dist']>mid_mag_out_radius)
        #print('far_pairs',far_pairs)
        #print("other_results[far_pairs]['dist']",other_results[far_pairs]['dist'])
        print('checking far pairs')
        far_bad_inds=np.where(mag_diff[far_pairs]>far_mag_func(other_results[far_pairs]['dist']))
        print('far pairs checked')
        print('far_bad_inds',far_bad_inds)
        #plt.scatter(other_results[far_pairs]['dist'],far_mag_func(other_results[far_pairs]['dist']),label='far boundary')
        if far_bad_inds[0].shape[0]>0:
            print('bad far stars')
            good_field=False
    #plt.scatter(other_results['dist'],mag_diff,label='field stars')
    #plt.legend()
    #plt.xlabel('distance (arcseconds)')
    #plt.ylabel('mag diff')
    #plt.show()
    if good_field:
        print('No magnitude issue stars in field based on separation-brightness functions')
    return good_field

#def get_panstarrs(source_id):
    #gaiadr2.panstarrs1_best_neighbour

    
    
    #return ps_data
#ps_table=Gaia.load_table('gaiadr2.panstarrs1_best_neighbour')
#print(ps_table)
#sys.exit()
print('\n\nMax number of stars in ',bright_search_radius, '" radius:', find_star_limit(bright_search_radius))
print('\n\nMax number of stars in ',sepbright_radius, '" radius:', find_star_limit(sepbright_radius))
dense_number=find_star_limit(bright_search_radius)

target_num = 1
for ra,dec,name in zip(ra_array, dec_array,name_array):
    begin_time = time.time()
    print("\n\ntarget ", str(target_num) + '/' + str(num_targets),'\n\n')
    results= cone_search(ra,dec)
    results['dist']=results['dist']*3600. #converting the distances into arcseconds
    print('\n\ndist in arcseconds:')
    print(results['dist'])
    #print(results['dist']*3600.)
    #print(results.shape)
    #new_col= Column(name, name='name')
    #print(results[0])
    #print(results[0])
    #output_row= results[0]
    name=str(name)
    #name='SDSS'+name
    print(name)
    print(type(name))
    good_field=True #this will be the boolean used to determine if we keep the results after this.
    try:
        table_length = len(results['dist'])
        print(table_length)
        if (filter_confused_sources and table_length>1):
            print("Too many sources in aperture")
            good_field=False
        #the rest of the checks only run if prior checks haven't ascertained that the field is bad.
        else:
            #print('checking other things')
            if ((filter_bright_contam or filter_dense) and (good_field)):
                bright_results= cone_search(ra,dec, search_radius=bright_search_radius)
                if (filter_bright_contam and np.nanmin(bright_results['phot_g_mean_mag']) < bright_star_limit):
                    print("Star brighter than G=" + str(bright_star_limit) + " , so likely contamination, thus excluding.")
                    good_field=False
                elif (filter_dense and len(bright_results['dist']) >= dense_number):
                    print("Stellar density greater than 100,000 deg^-2, there are ", len(bright_results['dist']), " stars")
                    good_field=False
                else:
                    pass
            if (filter_sepbright_func and good_field):
                print('\n\nGetting ',sepbright_radius,'arc second cone search')
                big_results=cone_search(ra,dec,search_radius=sepbright_radius)
                print('Cone Search completed.\n\n')
                big_results['dist']=big_results['dist']*3600.
                print('dist in arcseconds')
                print(big_results['dist'])
                good_field=check_distance_mag(big_results)
        if good_field:
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
                print('adding results to output table because field is good apparently...')
                collected_results.append(results[0])
                #collected_results.append(output_row)
                #collected_results.append(results)
            except IndexError:
                print("index error")
                pass
            
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
    except TypeError as error:
        print(error)
        pass

    else:
        pass
    end_time= time.time()
    print("\n\ntarget ", str(target_num) + '/' + str(num_targets),'\n\n')
    print('Target search time: ', find_time_difference(begin_time, end_time))
    print('\nTotal elapsed time', find_time_difference(start, end_time))
    print('\n\n\n===============')
    target_num += 1
    
    
stacked_results= vstack(collected_results)
stacked_results.pprint()

added_colnames=['pairdistance','sep_AU','R_chance_align','ra_wd','dec_wd','phot_g_mean_mag_wd']

stacked_results.add_columns(all_table[added_colnames])
stacked_results.pprint()

stacked_results.write(output_file, format='ascii.csv', overwrite=True)

end_overall= time.time()
print('\n\n')
print('Total elapsed time: ', find_time_difference(start, end_overall))

#for thing in coord_list:
    #print(thing.ra.to(u.hourangle), thing.dec)
