import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky

import pyregion
import pyregion._region_filter as filter
import pyregion._region_filter as region_filter
from astropy.wcs import WCS
from astropy.io.fits import Header

w = WCS("""
NAXIS   =                    2 / number of data axes
NAXIS1  =                 1000 / length of data axis
NAXIS2  =                 1000 / length of data axis
EQUINOX =                 J2000 / default
CTYPE1  = 'RA---TAN'
CRVAL1  =  1000
CRPIX1  =  1000
CDELT1  = -0.001
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--TAN'
CRVAL2  =  1000
CRPIX2  =  1000
CDELT2  =  0.001
CUNIT2  = 'deg     '
""")


loc = '/Users/Steven/Documents/prepare_simsack/'
mat = np.loadtxt(loc+'irac_all.cat', usecols=[0,1,8]) 
ra_IRAC = mat[:,0]
dec_IRAC = mat[:,1]
mag_IRAC = mat[:,2]
use = mag_IRAC <= 25.5
ra_IRAC = ra_IRAC[use]
dec_IRAC = dec_IRAC[use]
mag_IRAC = mag_IRAC[use]

c = SkyCoord(ra=ra_IRAC*u.degree, dec=dec_IRAC*u.degree) 
table = match_coordinates_sky(c, c, nthneighbor=2) 
table2 = np.array(table)
dist_arc = table2[1,:]*3600

SOURCE_LIST = np.append(np.linspace(10,18,21),np.linspace(18.2,27.0,45))

ra_rand =  np.min(ra_IRAC) + np.random.rand(100000)*(np.max(ra_IRAC)-np.min(ra_IRAC))
dec_rand =  np.min(dec_IRAC) + np.random.rand(100000)*(np.max(dec_IRAC)-np.min(dec_IRAC))


loc = '/Users/Steven/Documents/region_files/COSMOS_REGION/region-files/'
r = pyregion.open(loc+'cosmos_cen.reg').as_imagecoord(w)
in_range = np.zeros(np.size(ra_rand), dtype=bool)
r2 = r.get_filter(w)

for i in range(0,np.size(r)):
  print i
  r3 = r2[i]
  in_mask = r3.inside(x,y)
  in_range[in_mask] = True
print 'stop'
in_mask = in_range

print in_mask
