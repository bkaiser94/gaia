"""
Created by Ben Kaiser (UNC-Chapel Hill) 2023-06-06.

This should call custom_photometric_filters.py and the GaiaXPy function to make the calibrated
spectra in addition to whatever large file of targets with continuous XP spectra. It should then
produce a corresponding table of fluxes in each of the custom bands and then those should
in turn be used to produce informal color indices. I think I'll save the overall flux per band to the
general file and then I'll save a selection of a few of the color indices that can be created from 
the bands. I'll probably take the difference of the log of the flux in each band or something like 
that.




"""
