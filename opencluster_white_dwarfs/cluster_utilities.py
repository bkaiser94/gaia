"""
Created by Ben Kaiser (UNC-Chapel Hill) 2025-09-19

Help from Co-Pilot

This will contain some functions that I'll probably need to use more than once across scripts but that don't fit nicely into one of them.

The function that jumps to mind is selecting the relevant spectral types within the dataset.




"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from glob import glob
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column, vstack, join, MaskedColumn,hstack
#import scipy.interpolate as scinterp
import time
#import pyvo
#from astroquery.vizier import Vizier

#from astroquery.gaia import Gaia
#from astroquery.xmatch import XMatch


sys.path.append('../')

def limit_spectype(input_table,spectype='DA',spectype_colname='mwdd_spectype',flexible_spectype=True):
    
    """
    flexible_spectype: Boolean that dictates if the spectral type string must be an exact match or just contain the spectral type. E.g. "DA" with True will return "DA", "DAH", "DA:", "DAB", while False will only return "DA"
    
    
    """
    col = input_table[spectype_colname]

    if hasattr(col, 'mask'):
        unmasked = ~col.mask
    else:
        unmasked = np.ones(len(col), dtype=bool)
    subcol=col[np.where(unmasked)]
    subtable=input_table[np.where(unmasked)]
    #arr = col.astype(str)
    if flexible_spectype:

        
        #mask = np.char.find(col.astype(str),spectype) >= 0
        matches = np.char.find(np.array(subcol.astype(str)), spectype) #I just needed to convert to an array instead of running this on a column.
        #mask = unmasked & mask
        indices = np.where(matches>=0)
    else:
        indices=np.where(subcol==spectype)

    #output_table=input_table[indices].copy()
    output_table=subtable[indices].copy()
    print('Matched spectral types for '+spectype,', flexible_spectype',flexible_spectype,', count:',len(output_table))
    return output_table
