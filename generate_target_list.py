"""
Created by Ben Kaiser (UNC-Chapel Hill) 2019-01-09

Should take one of the fuller Gaia tables (meaning like all of the columns for the small target numbers), and
output a SOAR-compatible target list (so space-delimited or tab)

"""

from __future__ import print_function
import numpy as np
#from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table, vstack, Column

#input_file= '20190107_chris_merge_gaia.csv'
#input_file='20190405_purple_search_gaia_unique.csv'
#input_file='exc1_8_2_2_purple_search_gmaglimit_gaia_sc.csv'
#input_file='Lindegren_odd_survivors_gaia_sc.csv'
#input_file='sdssj1330_similar_subset_gaia_sc.csv'
#output_file= '20190107_chris_merge_targlist.txt'
#input_file='expanded_purple_search_gmaglimit_gaia_sc.csv'
#input_file= '20190516B_retargeted_purple_search_gaia_scbd.csv'
#input_file= '20190820B_priorityupdate_retargeted_purple_search_gaia_scbd.csv'
#input_file= '20190829B_priorityupdate_retargeted_purple_search_gaia_scbd.csv'
#input_file= '20191213B_priorityupdate_retargeted_purple_search_gaia_scbd.csv'
#input_file='20200110B_priorityupdate_retargeted_purple_search_gaia_scbd.csv'
#input_file='20200113B_priorityupdate_retargeted_purple_search_gaia_scbd.csv'
#input_file='20200712B_priorityupdate_retargeted_purple_search_gaia_scbd.csv'
#input_file= '20190516B_retargeted_purple_search_gaia_scbd_20210117_update.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20210707_update_neededobs.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20211117_update.csv'
#input_file='20190516B_retargeted_purple_search_gaia_scbd_20211205_update.csv'
#input_file='dimWDMS_F0toK7_eDR3.fits'
#input_file='20190917_alkaliWD_attempt2_gaia_scbd.csv'
#input_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf.csv'
#input_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_loosecutsonly_justWDradec_gaiaDR3_d_sbf.csv'
#input_file='WDs_mistaken_for_Kstars_nofehlimit_by_SDSS_gaiaDR3_wnames_andlabels.csv'
#input_file= 'josh_object.csv'
#input_file='dimWDMS_allMS_minsepfunc_eDR3_highconf_notrust_justWD_gaiaDR3_d_sbf_Vincent2023bClass_20241216update.csv'
#input_file='20250313_1125_hbetadipper_locus_500pcG18_gaiaDR3_20250313namesaddedfurther.csv'
#input_file='interesting_wds/gf21_GaiaeDR3_faintWDs_gaiaadded_veryinterestingwds_simbadadded.csv'
#input_file='interesting_wds/gf21_GaiaeDR3_faintWDs_gaiaadded_prettyinterestingwds_simbadadded.csv'
#input_file='Miller_clusterWDs_susages_gaiaDR3_.csv'
input_file='HR24members_GF21maincat_crossmatch_simbadadded_mwddadded_wdagesadded_agediffsig_gaiaadded_kaisernames_sigonly.csv'

#input_file='20190516B_retargeted_purple_subset.csv'
#comment_string='retarg_purple'
comment_string='clusterWD_sus_age_MWDDparams'
num=''
start_target_num=865
output_file= input_file.split('.')[0]+'_targlist.txt'
#output_file= input_file.split('.')[0]+'_400M1needed_targlist.txt'
#input_table= Table.read(input_file, format= 'ascii.csv')
input_table=Table.read(input_file)
delimiter='\t'

name_base='WDJ'
#name_base='GaiaJ'

sort_by_ra= True
list_length=input_table['ra'+num].shape[0]
target_num_array=np.arange(start_target_num,list_length+start_target_num+1,1)


input_table.pprint()
if sort_by_ra:
    sorted_order= np.argsort(input_table['ra'+num])
    print(sorted_order)
    print('sorted_order:')
    input_table=input_table[sorted_order]
    input_table.pprint() 

try:
    name_array= np.array(input_table['kcname'])
    #name_array= np.array(input_table['WDJname'])
except KeyError as error:
    print(error)
    print('No names column, so just making all Gaia names')
    name_array=np.full(input_table['ra'+num].shape, '.')

    #output_array= output_array[:,sorted_order]

coords= coord.SkyCoord(input_table['ra'+num], input_table['dec'+num], unit=(u.deg, u.deg))
string_coords= coords.to_string(style='hmsdms')
replace_chars= ['d','h','m']
ra_list=[]
dec_list=[]
epoch_list=[]
mag_list=[]
name_list=[]
mag_string='g'
full_mag_string= 'phot_'+mag_string+'_mean_mag'+num
#full_mag_string='Gmag'
#for thing, name, mag  in zip(string_coords, name_array, input_table[full_mag_string]):
#for thing, name, mag, target_num  in zip(string_coords, name_array, input_table[full_mag_string], input_table['target_num']):


    

for thing, name, mag, target_num  in zip(string_coords, name_array, input_table[full_mag_string], target_num_array):
#for thing, name, mag, target_num  in zip(string_coords, name_array, input_table['app_sp_type'], target_num_array):


    print(name)
    print(thing)
    for char in replace_chars:
        thing= thing.replace(char, ':')
    thing=thing.replace('s', '')
    split_string=thing.split(' ')
    ra= split_string[0][:11] #limiting precision of the decimal
    dec= split_string[1][:12] #limiting precision of the decimal, need the additional index for sign
    small_ra= ra.replace(':', '')[:4]
    small_dec=dec.replace(':','')[:5] #need additional index for + or -
    if name=='':
        name=name_base+small_ra+small_dec
        print('new name:',name)
    if name=='--':
        name=name_base+small_ra+small_dec
        print('new name:',name)
    elif name=='.':
        #name='GaiaJ'+small_ra+small_dec
        name=name_base+small_ra+small_dec
        print('new name:', name)
    elif name=='none':
        #name='GaiaJ'+small_ra+small_dec
        name=name_base+small_ra+small_dec
        print('new name:', name)
    elif name == '0':
        #name='GaiaJ'+small_ra+small_dec
        name=name_base+small_ra+small_dec
        print('new name:', name)
    print('thing',thing)
    #remove spaces from names
    name=name.replace(' ','')
    name=str(int(target_num))+'_' + name
    
    name_list.append(name)
    #ra=ra.replace(':',' ')
    #dec=dec.replace(':',' ')
    
    
    ra_list.append(ra)
    dec_list.append(dec)
    epoch_list.append('2000.0')
    #mag_list.append(mag_string.upper()+'='+str(mag)[:4])
    mag_list.append(mag_string.upper()+'='+str(mag)[:4]+','+comment_string)
    #mag_list.append(mag)
    #for other in thing:
        #print(other.replace(['d','h','m','s'], ':'))
print(name_array)

output_array= np.vstack([name_list, ra_list, dec_list, epoch_list, mag_list]).T
#priority_good= np.where(input_table['priority']< 10)
#priority_good= np.where(input_table['400m2_need_bool']== 1)
#priority_good= np.where((input_table['400m1_need_bool']!= 0) and  (input_table['400m2_need_bool']!=0))
#priority_good= np.where(input_table['400m1_need_bool']== 1)
#priority_good= np.where((input_table['400m1_need_bool']== 1) or (input_table['400m2_need_bool']==1))
#priority_good= np.where((input_table['obs_need_bool']== 1))
#output_array=output_array[priority_good]
output_array=output_array.T
#if sort_by_ra:
    #sorted_order= np.argsort(output_array[1])
    #print(sorted_order)
    #print('sorted_order:',output_array[:,sorted_order])
    #output_array= output_array[:,sorted_order]
print(output_array)
print(output_array.dtype)
output_array= output_array.T
np.savetxt(output_file,output_array,  delimiter= delimiter, fmt= '%1s')


