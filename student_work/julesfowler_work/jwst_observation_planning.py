## -- OBSERVATION PLANNING WITH JWST 
""" This is pretty simple, contains the code to read in the Sonora model,
convert the flux, write out the new file, and print out the normalization we
need for APT."""

## -- IMPORTS
from astropy.io import ascii
import numpy as np

## -- GO
# Read in file 
data = np.loadtxt('../../Project 5--JWST/sp_t350g100nc_m0.0', skiprows=2)
wv, flux = data[:, 0], data[:, 1]

# Calculate flux conversion factor
r_jup = 69.911e6
pc = 3.086e16

r_wise = 1*r_jup
d_wise = 5.7*pc

flux_convert = (r_wise**2)/(d_wise**2)
flux = flux*flux_convert*1e23

# Trim the file to 2-6 um so it's not too big
flux = flux[wv > 2]
wv = wv[wv > 2]
flux = flux[wv < 6]
wv = wv[wv < 6]

# Make a new dictionary
data_out = {'#wavelength': wv, '#flux': flux}
ascii.write(data_out, 'wise_out.dat')

# Find peak to renormalize to
index = np.where(np.max(flux) == flux)[0][0]
print(f'Peak flux is {flux[index]} at lambda={wv[index]}')



