"""
Created by Ben Kaiser (UNC- Chapel Hill) 2018-06-08

Take several input parameters (such as being located within 100 pc), and query the Gaia database for DR2 data.
Then return the output to a text file that can be read by other programs, such as one to make a color-magnitude
diagram.

"""
import numpy as np
import astroquery

