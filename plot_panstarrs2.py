"""
Created by Ben Kaiser (UNC-Chapel Hill) 2022-05-02.

Hopefully this will plot the panstarrs data for the purple objects that I retrieved from MAST.

Probably is going to have to be run in python3

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
#import astropy


#import passband_model_convolution as pmc
import gaia_extinction
#import wdatmos
import plotting_dicts as pod


#color_string_base='MeanPSFMag'
color_string_list=[
    'g',
    'r',
    'i',
    'z',
    'y'
    
    ]
    
color_string_base='PSFMag'
#color_string_list=[
    #'r',
    #'i',
    #'z'
    #]

#input_file='20200819_retargeted_purple_search_wpanstarrs_REALLYjust_ra_dec_csv_5_2_2022_panstarrs2_params_sortbyra_goodmags.csv'
#input_file='20200819_retargeted_purple_search_wpanstarrs_REALLYjust_ra_dec_csv_5_2_2022_panstarrs2_params.csv'
#input_file='20200819_retargeted_purple_search_wpanstarrs_REALLYjust_ra_dec_csv_5_3_2022_panstarrs1_params.csv'
input_file='20200819_retargeted_purple_search_wpanstarrs_REALLYjust_ra_dec_csv_5_3_2022_panstarrs1stack_params.csv'

input_table=Table.read(input_file)



color_lists=[
    
    
    
    
    
    
    ]

####################

def clean_table(input_table):
    output_table=input_table.copy()
    print('\n\nCleaning table.')
    print('Starting from', len(output_table))
    for color in color_string_list:
        print('Cleaning ', color)
        #good_inds=np.where(output_table[color+color_string_base]>-990.)
        #output_table=output_table[good_inds]
        #good_inds=np.where(output_table[color+color_string_base+'Err']>-990.)
        #output_table=output_table[good_inds]
        
        print('New count:', len(output_table))
    
    sorted_order=np.argsort(output_table['_ra_'])
    output_table=output_table[sorted_order]
    
    
    
    return output_table


################3

input_table=clean_table(input_table)

output_file_parts=input_file.split('.')[0]
output_file=output_file_parts+'_properRAsort_goodmags.csv'
input_table.write(output_file, format='ascii.csv')


def plot_panstarrs_colors(c1c2='g-r', c3c4='i-z', input_table=input_table,color='',label=''):
    c1,c2=c1c2.split('-')
    c3,c4=c3c4.split('-')
    if color=='':
        plt.scatter(input_table[c1+color_string_base]-input_table[c2+color_string_base], input_table[c3+color_string_base]-input_table[c4+color_string_base],label=label)
    else:
        plt.scatter(input_table[c1+color_string_base]-input_table[c2+color_string_base], input_table[c3+color_string_base]-input_table[c4+color_string_base], color=color,label=label)
    plt.xlabel(c1c2)
    plt.ylabel(c3c4)
    plt.show()
    return


def plot_colors(c1c2='g-r', c3c4='i-z', input_table=input_table,color='',label=''):
    """
    The reformatted strings that I used in the gaia files should be used here
    
    """
    c1,c2=c1c2.split('-')
    c3,c4=c3c4.split('-')
    #ps_string='ps2_'
    ps_string='ps1_'
    if color=='':
        plt.scatter(input_table[ps_string+c1+'_mean_mag']-input_table[ps_string+c2+'_mean_mag'], input_table[ps_string+c3+'_mean_mag']-input_table[ps_string+c4+'_mean_mag'],label=label)
    else:
        plt.scatter(input_table[ps_string+c1+'_mean_mag']-input_table[ps_string+c2+'_mean_mag'], input_table[ps_string+c3+'_mean_mag']-input_table[ps_string+c4+'_mean_mag'], color=color,label=label)
    plt.xlabel(c1c2)
    plt.ylabel(c3c4)
    plt.show()
    return

def plot_colors_vs_ra(c1c2='g-r'):
    c1,c2=c1c2.split('-')

    plt.scatter(input_table['MatchRA'], input_table[c1+color_string_base]-input_table[c2+color_string_base])
    plt.xlabel('Matched RA')
    plt.ylabel(c1c2)
    plt.show()
    return

if __name__ =='__main__':
    plot_panstarrs_colors(c1c2='g-r', c3c4='i-z')
    plot_panstarrs_colors(c1c2='r-i', c3c4='i-z')

    plot_panstarrs_colors(c1c2='g-r',c3c4='r-i')
    plot_panstarrs_colors(c1c2='r-i',c3c4='i-y')


    plot_colors_vs_ra(c1c2='g-r')
    plot_colors_vs_ra(c1c2='r-i')
    plot_colors_vs_ra(c1c2='i-z')
    plot_colors_vs_ra(c1c2='z-y')









