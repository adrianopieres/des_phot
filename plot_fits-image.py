# -*- coding: utf-8 -*-

import numpy as np
import astropy.io.fits as fits
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
cmap = mpl.cm.get_cmap("gray")
from matplotlib import colors
cmap.set_bad('k')

x_min, x_max, y_min, y_max = 0, 1000, 0, 1000

image_data = fits.getdata('DES0221-0958_g.fits', ext=0)
f = fits.open('DES0221-0958_g.fits')
w = WCS(f[1].header)

table = fits.getdata("cats/DES0221-0958_gold_y6.fits")
RA = table['ra']
DEC = table['dec']

sky = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)

x, y = w.world_to_pixel(sky)

im = plt.imshow(image_data[x_min:x_max, y_min:y_max], cmap='gray', norm=colors.LogNorm(vmin=1., vmax=np.max(image_data)), origin='lower')

plt.scatter(x, y, color='None', edgecolor='b')

plt.xlim([x_min, x_max])
plt.ylim([y_min, y_max])
plt.colorbar(im, cmap=cmap)
plt.show()
