"""
Created by Ben Kaiser (UNC-Chapel Hill) 2018-07-11

This should calculate the extinction in each Gaia passband by first calculating the reddening of the object and
then using some method (probably just the standard one from DR2HRD) to get the A_0 extinction to then compute the extinction in each passband.

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
import convert_colours

#extinction coefficients from DR2HRD 2018 Table 1
#the 0th entry is blank to match the indexing of the table
c_coeffs= {'g': [0, 0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0094],
           'bp': [0, 1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043],
           'rp': [0, 0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005,0.0006]}

def get_reddening(obs_bp_rp, model_bp_rp):
    obs_b_v = convert_colours.find_B_V(obs_bp_rp)
    model_b_v= convert_colours.find_B_V(model_bp_rp)
    return model_b_v- obs_b_v

def get_a_0(reddening):
    return 3.1*reddening

def get_a_x(obs_bp_rp,  model_bp_rp,passband_string = 'g'):
    """
    Equation 1 from DR2HRD 2018
    """
    reddening= get_reddening(obs_bp_rp, model_bp_rp)
    a_0= get_a_0(reddening)
    c_vals = c_coeffs[passband_string]
    k_x=1
    k_x= c_vals[1]+c_vals[2]*bp_rp+c_vals[3]*bp_rp**2+c_vals[4]*bp_rp**3+c_vals[5]*a_0+c_vals[6]*a_0**2+c[7]*bp_rp*a_0
    a_x = a_0*k_x
    return a_x


