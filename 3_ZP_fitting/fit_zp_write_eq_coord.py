import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from astropy.io import fits
from astropy.wcs import WCS
import glob
import matplotlib.pyplot as plt

def write_ZP(infile):
    """_summary_

    Parameters
    ----------
    infile : _type_
        _description_
    infile_DES : _type_
        _description_
    band : _type_
        _description_
    mag_lim_sat_DES : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    # Query example on easyaccess:
    # SELECT ra, dec, wavg_mag_psf_g, wavg_mag_psf_r, wavg_mag_psf_i, wavg_mag_psf_z, wavg_mag_psf_y, flags_gold, ext_coadd, TILENAME FROM Y6_GOLD_2_0 WHERE ABS(ext_coadd) < 2 AND flags_gold = 0 AND ABS(wavg_mag_psf_g) < 28. AND TILENAME = 'DES0221-0958';> DES0221-0958_gold_y6.fits
    image_DES = '../1_DES_images_and_cats/images/' + infile[16:-16] + 'fits.fz'
    infile_DES = '../1_DES_images_and_cats/cats/' + infile[16:-18] + 'gold_y6.fits'
    band = infile[29:30]
    f = fits.open(image_DES)
    w = WCS(f[1].header)
    hdu = fits.open(infile_DES, memmap=True)
    RA_DES = hdu[1].data.field('ra')
    DEC_DES = hdu[1].data.field('dec')
    MAG_DES = hdu[1].data.field('wavg_mag_psf_'+ band)
    hdu.close()
    
    cond = (MAG_DES > mag_lim_sat_DES[band] + 1.)&(MAG_DES < mag_lim_sat_DES[band] + 3.)
    
    RA_DES, DEC_DES, MAG_DES = RA_DES[cond], DEC_DES[cond], MAG_DES[cond]
    
    IDX, X_DAO, Y_DAO, MAG_DAO, MAG_ERR, SHARPNESS, CHI = np.loadtxt(infile, usecols=(0,1,2,3,4,5,6), unpack=True)
    
    sky = w.pixel_to_world(X_DAO, Y_DAO)
    RA_DAO = sky.ra.degree
    DEC_DAO = sky.dec.degree
    C_DAO = SkyCoord(ra=RA_DAO*u.degree, dec=DEC_DAO*u.degree)
    C_DES = SkyCoord(ra=RA_DES*u.degree, dec=DEC_DES*u.degree)

    idx_dao, idx_des, d2d, d3d = C_DES.search_around_sky(C_DAO, 1*u.arcsec)

    # plot Mag fotometria x Mag DES:
    # for i,j in zip(idx_dao, idx_des):
    #    print(C_DAO[i].ra.deg, C_DAO[i].dec.deg, C_DES[j].ra.deg, C_DES[j].dec.deg)

    best_fit = np.polyfit(
        MAG_DAO[idx_dao], MAG_DES[idx_des], deg=1) #, w=1/MAG_ERR[idx_dao]**2)

    ZP = (np.sum(MAG_DAO[idx_dao] - MAG_DES[idx_des]) / len(idx_des))

    xx = np.arange(10, 26, 0.5)
    plt.scatter(MAG_DAO[idx_dao], MAG_DES[idx_des], label='data', s=0.2)
    plt.plot(xx + ZP, xx, color='k', label='Fit')
    plt.xlabel('MAG in {} band'.format(band))
    plt.legend()
    plt.xlabel('Magnitude in daophot')
    plt.ylabel('Magnitude DES')
    plt.title('Best-Fitting for {} band: {:.4f}'.format(band,
              ZP))
    plt.legend(loc=2)
    plt.xlim([np.min(MAG_DAO[idx_dao]) - 0.5, np.max(MAG_DAO[idx_dao]) + .5])
    plt.ylim([mag_lim_sat_DES[band] +0.5, mag_lim_sat_DES[band] + 3.5])
    plt.savefig(infile[16:-18] + band +'_ZP' + '.png')
    plt.close()
    
    g = open(infile[16:-16] + '_calibrated_with_wcs.dat', 'w')
    print(infile[16:-16] + '_calibrated_with_wcs.dat')
    for i in range(len(IDX)):
        print(IDX[i], X_DAO[i], Y_DAO[i], RA_DAO[i], DEC_DAO[i], MAG_DAO[i] + ZP, MAG_ERR[i], SHARPNESS[i], CHI[i], file=g)
    g.close()

pre_cal_cat = glob.glob('../2_photometry/*pre_cal.dat')

mag_lim_sat_DES = {'g': 16.2,
		   'r': 16.7,
		   'i': 16.8,
		   'z': 15.5,
		   'Y': 15.6}


bands = ['g', 'r', 'i', 'z', 'Y']

for i in pre_cal_cat:
    print(i)
    write_ZP(i)

