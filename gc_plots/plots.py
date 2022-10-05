import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from itertools import compress
from matplotlib.patches import Polygon

ra_des_ftp, dec_des_ftp = np.loadtxt('des-round17-poly.txt', unpack=True)
ra_des_ftp[ra_des_ftp > 180.] -= 360.

L, B = np.loadtxt('Harris.dat', usecols=(1,2), unpack=True)
GC_name = np.loadtxt('Harris.dat', usecols=(0), dtype=str, unpack=True)

c_gal = SkyCoord(l=L*u.degree, b=B*u.degree, frame='galactic')

ra_gc, dec_gc = c_gal.icrs.ra.deg, c_gal.icrs.dec.deg

ra_gc[ra_gc > 180.] -= 360.

ra_tile_c, dec_tile_c, rall, decll, raul, decul, raur, decur, ralr, declr = np.loadtxt('coaddtiles-20121015.csv', delimiter=',', usecols=(3,4,9,10,11,12,13,14,15,16), unpack=True)

ra_tile_c[ra_tile_c > 180.] -= 360.

tile_name = np.loadtxt('coaddtiles-20121015.csv', usecols=(2), delimiter=',', dtype=str, unpack=True)

polygon = Polygon([(i,j) for i,j in zip(ra_des_ftp, dec_des_ftp)])
inside = [polygon.contains(Point(i, j))[0] for i,j in zip(ra_tile_c, dec_tile_c)]

ra_tile_c[ra_tile_c > 180.] -= 360.
rall[rall > 180.] -= 360.
raul[raul > 180.] -= 360.
raur[raur > 180.] -= 360.
ralr[ralr > 180.] -= 360.

ra_tile_c = list(compress(ra_tile_c, inside))
dec_tile_c = list(compress(dec_tile_c, inside))
rall = list(compress(rall, inside))
decll = list(compress(decll, inside))
raul = list(compress(raul, inside))
decul = list(compress(decul, inside))
raur = list(compress(raur, inside))
decur = list(compress(decur, inside))
ralr = list(compress(ralr, inside))
declr = list(compress(declr, inside))
tile_name = list(compress(tile_name, inside))
fig,ax = plt.subplots()

for i in range(len(rall)):
    y = np.array([[rall[i], decll[i]], [raul[i], decul[i]], [raur[i],decur[i]], [ralr[i],declr[i]], [rall[i], decll[i]]])
    p = Polygon(y, edgecolor = 'k', alpha=0.2)
    ax.add_patch(p)
    # inside2 = [p.contains(Point(ii, jj)) for ii,jj in zip(ra_gc, dec_gc)]
    # if any(inside2):
    #     plt.text(ra_tile_c[i], dec_tile_c[i], tile_name[i], fontsize=6, color='grey')

for i in range(len(rall)):
    y = np.array([[rall[i], decll[i]], [raul[i], decul[i]], [raur[i],decur[i]], [ralr[i],declr[i]], [rall[i], decll[i]]])
    p = Polygon(y, edgecolor = 'k', alpha=0.2)
    inside2 = [p.contains(Point(ii, jj))[0] for ii,jj in zip(ra_gc, dec_gc)]
    print(inside2)
    if any(inside2):
        plt.text(ra_tile_c[i], dec_tile_c[i], tile_name[i], fontsize=12, color='k')

for i in range(len(ra_gc)):
    plt.text(ra_gc[i], dec_gc[i], GC_name[i], fontsize=12)

plt.plot(ra_des_ftp, dec_des_ftp, color='r')
plt.scatter(ra_gc, dec_gc, marker='o', c='k', s=0.2)
ax.set_xlim([100,-65])
ax.set_ylim([-68, 7])
plt.savefig('ftp.png')
plt.close()
