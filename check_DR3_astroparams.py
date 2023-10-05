"""
Created by Ben Kaiser (UNC-Chapel Hill) 2023-09-28


Query for astrophysical parameters from the MCMC Multiple Star Classifier in Gaia DR3.

Will probably use the source_id to do this matching since that's what appears to be present in both tables.




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


input_file='WDJ1948m1011_gaiaDR3.csv'
