import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
import glob
from astropy.io.fits import getdata
import itertools

tiles = glob.glob('../3_ZP_fitting/DES*_g._calibrated_with_wcs.dat')
bands = ['g', 'r', 'i', 'z', 'Y']

for i in tiles:
    idx_all = []
    
    id_g, x_g, y_g, ra_g, dec_g, mag_g, magerr_g, sharp_g, chi_g = np.loadtxt(i.replace('_g', '_g'), unpack=True)
    id_g = [int(i) for i in id_g]
    id_r, x_r, y_r, ra_r, dec_r, mag_r, magerr_r, sharp_r, chi_r = np.loadtxt(i.replace('_g', '_r'), unpack=True)
    id_r = [int(i) for i in id_r]
    id_i, x_i, y_i, ra_i, dec_i, mag_i, magerr_i, sharp_i, chi_i = np.loadtxt(i.replace('_g', '_i'), unpack=True)
    id_i = [int(i) for i in id_i]
    id_z, x_z, y_z, ra_z, dec_z, mag_z, magerr_z, sharp_z, chi_z = np.loadtxt(i.replace('_g', '_z'), unpack=True)
    id_z = [int(i) for i in id_z]
    id_Y, x_Y, y_Y, ra_Y, dec_Y, mag_Y, magerr_Y, sharp_Y, chi_Y = np.loadtxt(i.replace('_g', '_Y'), unpack=True)
    id_Y = [int(i) for i in id_Y]

    idx_ = list(itertools.chain(id_g, id_r, id_i, id_z, id_Y))
    idx_un = sorted(set(idx_))
    
    with open(i[16:].replace('_g._calibrated_with_wcs.dat', '_grizY_calibrated.dat'), 'w') as outfile:
        for i in idx_un:
            # print(i, file=outfile)

            if i in id_g:
                arg = id_g.index(i)
                print('{:d} {:.3f} {:.3f} {:.5f} {:.5f} {:.5f} {:.5f} {:.3f} {:.3f}'.format(id_g[arg], x_g[arg], y_g[arg], ra_g[arg], dec_g[arg], mag_g[arg], magerr_g[arg], sharp_g[arg], chi_g[arg]), file=outfile, end=' ')
            else:
                print('-999 -999.999 -999.999 -999.999 -999.999 -999.99999 -999.999 -999.999 -999.999', file=outfile, end=' ')


            if i in id_r:
                arg = id_r.index(i)
                print('{:d} {:.3f} {:.3f} {:.5f} {:.5f} {:.5f} {:.5f} {:.3f} {:.3f}'.format(id_r[arg], x_r[arg], y_r[arg], ra_r[arg], dec_r[arg], mag_r[arg], magerr_r[arg], sharp_r[arg], chi_r[arg]), file=outfile, end=' ')
            else:
                print('-999 -999.999 -999.999 -999.999 -999.999 -999.99999 -999.999 -999.999 -999.999', file=outfile, end=' ')
    

            if i in id_i:
                arg = id_i.index(i)
                print('{:d} {:.3f} {:.3f} {:.5f} {:.5f} {:.5f} {:.5f} {:.3f} {:.3f}'.format(id_i[arg], x_i[arg], y_i[arg], ra_i[arg], dec_i[arg], mag_i[arg], magerr_i[arg], sharp_i[arg], chi_i[arg]), file=outfile, end=' ')
            else:
                print('-999 -999.999 -999.999 -999.999 -999.999 -999.99999 -999.999 -999.999 -999.999', file=outfile, end=' ')


            if i in id_z:
                arg = id_z.index(i)
                print('{:d} {:.3f} {:.3f} {:.5f} {:.5f} {:.5f} {:.5f} {:.3f} {:.3f}'.format(id_z[arg], x_z[arg], y_z[arg], ra_z[arg], dec_z[arg], mag_z[arg], magerr_z[arg], sharp_z[arg], chi_z[arg]), file=outfile, end=' ')
            else:
                print('-999 -999.999 -999.999 -999.999 -999.999 -999.99999 -999.999 -999.999 -999.999', file=outfile, end=' ')


            if i in id_Y:
                arg = id_Y.index(i)
                print('{:d} {:.3f} {:.3f} {:.5f} {:.5f} {:.5f} {:.5f} {:.3f} {:.3f}'.format(id_Y[arg], x_Y[arg], y_Y[arg], ra_Y[arg], dec_Y[arg], mag_Y[arg], magerr_Y[arg], sharp_Y[arg], chi_Y[arg]), file=outfile)
            else:
                print('-999 -999.999 -999.999 -999.999 -999.999 -999.99999 -999.999 -999.999 -999.999', file=outfile)
    outfile.close()
