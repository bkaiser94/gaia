"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-07-11

This should perform the Gaia colour conversions outlined in Appendix A of Evans et al. 2018 (Photometric content and validation). 

The coefficients are taken from Table A.2 of Evans et al. 2018, and the bounds of the fit are contained in Table A.1.
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



#import wdatmos
import spec_plot_tools as spt


colour_coeffs= {
    'Johnson-Cousins':{
        "G-V":{
            "B-V": {"coeffs":[-0.02907,-0.02385,-0.2297,-0.001768], 
                    "bounds":[-0.3, 2.4],
                    "sigma":0.06285},
            "BP-RP":{"coeffs":[-0.01760,-0.006860,-0.1732,0],
                     "bounds":[-0.5,2.75],
                     "sigma":0.045858}}},
        "Hipparcos":{
            "G-Hp":{
                "B-V":{"coeffs":[-0.03704,-0.3915,0.01855,-0.03239,],
                       "bounds":[-0.2,1.5],
                       "sigma":0.02502},
                "BP-RP":{"coeffs":[-0.01968,-0.2344,-0.1200,0.01490],
                         "bounds":[-0.5,4.0],
                         "sigma":0.06875}}}}
                
def calc_colour(input_colour, coeffs):
    return coeffs[0] + coeffs[1]*input_colour+coeffs[2]*input_colour**2 + coeffs[3]*input_colour**3


def find_B_V(input_colour, filter_set = 'Johnson-Cousins',start_colour= 'BP-RP',end_colour = 'B-V'):
    """
    Can't handle arrays as inputs....
    """
    if filter_set == 'Johnson-Cousins':
        mid_colour = 'G-V'
    elif filter_set == 'Hipparcos':
        mid_colour = 'G-Hp'
    start_dict= colour_coeffs[filter_set][mid_colour][start_colour]
    start_constant = calc_colour(input_colour, start_dict['coeffs'])
    if ((input_colour<=start_dict['bounds'][0]) or (input_colour >= start_dict['bounds'][1])):
            print(start_colour + ':' + str(input_colour) + ' outside bounds: ' + str(start_dict['bounds']) + ' for table')
    end_dict = colour_coeffs[filter_set][mid_colour][end_colour]
    def func_to_solve(output_colour, coeffs = end_dict['coeffs'], start_constant= start_constant):
        return coeffs[0] + coeffs[1]*output_colour+coeffs[2]*output_colour**2 + coeffs[3]*output_colour**3-start_constant
    soln_object = sciop.root(func_to_solve, 0)
    soln = soln_object.x
    if ((soln<=end_dict['bounds'][0]) or (soln >= end_dict['bounds'][1])):
            print(end_colour + ':' + str(soln) + ' outside bounds: ' + str(end_dict['bounds']) + ' for table')
    print("soln:", soln)
    return soln
    
#find_B_V(0.688)
find_B_V(0.470)








