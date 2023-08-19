from __future__ import print_function
import subprocess
import numpy as np
from math import sqrt
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import pyraf
from pyraf import iraf
from pyraf.iraf import digiphot, daophot, imarith
import glob
import os
import matplotlib.pyplot as plt


def zp_from_file(tile):
    tiles = np.loadtxt('ZPs.dat', usecols=(0), dtype=str, unpack=True)
    zps = np.loadtxt('ZPs.dat', usecols=(2), unpack=True)
    idx = numpy.where(tiles == tile)
    return zps[idx]


def add_star(im_name, photfile, psfimage, addimage, minmag, maxmag, nstar):
    iraf.addstar(image=im_name + '[0]', photfile=photfile, psfimage=psfimage, addimage="default", minmag=minmag, maxmag=maxmag, nstar=nstar, verify='no', verbose='no')
    print('task addstar concluida em {}'.format(im_name))
    
    
def matching_stars(infile_phot, infile_DES):
    """_summary_

    Parameters
    ----------
    infile_phot : _type_
        _description_
    infile_DES : _type_
        _description_
    """
    RA, DEC, MAG_I, MAGI_ERR = np.loadtxt(infile_phot, usecols=(1, 2, 3, 4), unpack=True)
    RA_DES, DEC_DES, MAGG_DES, MAGG_ERR_DES, MAGR_DES, MAGR_ERR_DES, MAGI_DES, MAGI_ERR_DES, MAGZ_DES, MAGZ_DES_ERR, MAGY_DES, MAGY_ERR_DES = np.loadtxt(infile_DES, usecols=(2, 3, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17), delimiter=',', skiprows=1, unpack=True)

    # Agora configurando cada catalogo
    C_DAO = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
    C_DES = SkyCoord(ra=RA_DES*u.degree, dec=DEC_DES*u.degree)

    idx_dao, idx_des, d2d, d3d = C_DES.search_around_sky(C_DAO, 1*u.arcsec)

    # agora escrevendo ambas as coordenadas em ambos os arquivos:
    for i, j in zip(idx_dao, idx_des):
        print(C_DAO[i].ra.deg, C_DAO[i].dec.deg, MAG_I[i], MAGI_ERR[i], C_DES[j].ra.deg, C_DES[j].dec.deg, MAGI_DES[j], MAGI_ERR_DES[j], file=f)
    f.close()


def sky_sigma(im_data):
    """_summary_

    Parameters
    ----------
    im_data : array
        Data with counts about sky.

    Returns
    -------
    float, float
        Measurements about sky and standard deviation of sky.
    """
    sky = np.median(im_data)
    sigma = np.std(im_data[im_data < sky-np.min(im_data)])
    return sky, sigma


def detection(im_name):
    """Run script on the detections image.

    Parameters
    ----------
    im_name : str
        Data containing image name
    """
    im_data = fits.getdata(im_name, ext=0)

    sky, sigma = sky_sigma(im_data)

    iraf.daopars.unlearn
    iraf.photpars.unlearn
    iraf.datapars.unlearn
    iraf.fitskypars.unlearn
    iraf.centerpars.unlearn
    iraf.findpars.unlearn
    iraf.photpars.zmag = 25.00
    # Procurar no ds9 #Maximum equivalent exposure time (s)
    iraf.datapars.exposure = "EXPTIME"
    iraf.daopars.functio = "AUTO"
    iraf.datapars.gain = "GAIN"  # Maximum equivalent gain (e-/ADU)
    iraf.datapars.scale = 1
    # iraf.datapars.datamin = -10.0
    # iraf.datapars.datamax = 1500.0  # ou utilizar #1500
    iraf.daopars.saturated = 'yes'
    # de acordo com a imagem #desvio padra da #imexamin(comando 'r') no ds9
    iraf.datapars.fwhmpsf = 4.00
    iraf.datapars.sigma = sigma  # de acordo com a imagem
    iraf.datapars.readnoi = 0.0
    iraf.daopars.psfrad = 30.0  # de acordo com imagem # raio da maior imagem
    iraf.datapars.epadu = 1.  # de acordo com imagem
    iraf.daopars.fitrad = 6.  # de acordo com a imagem #(~1.5 ou 2.5 fwhm)
    iraf.daopars.recenter = 'yes'
    iraf.daopars.fitsky = 'yes'
    iraf.daopars.groupsk = 'yes'
    iraf.daopars.sannulu = 60.0
    iraf.fitskypars.annulus = 60
    iraf.fitskypars.dannulu = 10
    iraf.fitskypars.skyvalu = sky
    iraf.photpars.apertur = 12.0
    iraf.centerpars.cbox = 3
    iraf.centerpars.cthresh = 3.0

    iraf.daofind(image=im_name + '[0]', threshold=3.00,
                 output='default', verbose='no', verify='no')
    iraf.pdump(infile=im_name + '0.coo.1', fields="XCENTER, YCENTER, MAG, SHARPNESS, SROUND, GROUND",
               expr="(MAG != INDEF)&&(SHARPNESS != INDEF)&&(SROUND != INDEF)&&(GROUND != INDEF)", Stdout=im_name + '_det.dat')
    iraf.pdump(infile=im_name + '0.coo.1', fields="ID, XCENTER, YCENTER, MAG, SHARPNESS, SROUND, GROUND",
               expr="(MAG != INDEF)&&(SHARPNESS != INDEF)&&(SROUND != INDEF)&&(GROUND != INDEF)", Stdout=im_name + '_det.dat_with_n')


def ap_phot(im_name, coord_file):
    """ This function runs aperture photometry based on daophot phot task

    Parameters
    ----------
    im_name : str
        Name of the image.
    coord_file : str
        File name with the x and y coordinates of detected sources.
    max_n_psf : int
        Maximum number of psf candidates.
    """
    im_data = fits.getdata(im_name, ext=0)
    sky, sigma = sky_sigma(im_data)

    iraf.datapars.sigma = sigma
    iraf.fitskypars.skyvalu = sky

    iraf.phot(image=im_name+'[0]', coord=coord_file,
              output='default', verbose='no', verify='no')
    iraf.pdump(infile=im_name + '0.mag.1', fields="ID,XCENTER,YCENTER,MAG,MERR,SHARPNESS,CHI",
               expr="(MAG != INDEF)&&(MERR != INDEF)", Stdout=im_name + '0.mag')
    print('Ap_phot run for {}'.format(im_name))


def matching_stars(infile_phot, infile_DES):
    """ Matches the stars in infle and Dark Energy Survey.

    Parameters
    ----------
    infile_phot : str
        File name with the photometry.
    infile_DES : str
        File name with information from DES.
    """
    RA, DEC, MAG_I, MAGI_ERR = np.loadtxt(
        infile_phot, usecols=(1, 2, 3, 4), unpack=True)
    RA_DES, DEC_DES, MAGG_DES, MAGG_ERR_DES, MAGR_DES, MAGR_ERR_DES, MAGI_DES, MAGI_ERR_DES, MAGZ_DES, MAGZ_DES_ERR, MAGY_DES, MAGY_ERR_DES = np.loadtxt(
        infile_DES, usecols=(2, 3, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17), delimiter=',', skiprows=1, unpack=True)

    C_DAO = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
    C_DES = SkyCoord(ra=RA_DES*u.degree, dec=DEC_DES*u.degree)

    idx_dao, idx_des, d2d, d3d = C_DES.search_around_sky(C_DAO, 1*u.arcsec)

    for i, j in zip(idx_dao, idx_des):
        print(C_DAO[i].ra.deg, C_DAO[i].dec.deg, MAG_I[i], MAGI_ERR[i],
              C_DES[j].ra.deg, C_DES[j].dec.deg, MAGI_DES[j], MAGI_ERR_DES[j], file=f)
    f.close()


def min_dist(X, Y, MAG, RANGE_MAG, INIT_ANG_DIST=20.):
    """ Calculates the minimum distance in pixels.

    Parameters
    ----------
    X : array
        X cartesian coordinates of sources.
    Y : array
        Y cartesian coordinates of sources.
    MAG : array
        Array with the magnitudes of the stars.
    RANGE_MAG : float
        Range in magnitude to calculate minimum dstance.
    INIT_ANG_DIST : _type_, optional
        Minimum distance in pixels, by default 20.

    Returns
    -------
    array
        Minimum distance for stars. If the star is free from neighbours,
        distance is equal to 99.999
    """
    min_dist = np.repeat(99.999, len(X))
    for i in range(len(X)):
        cond = (X > X[i] - INIT_ANG_DIST) & (X < X[i] + INIT_ANG_DIST) & \
            (Y > Y[i] - INIT_ANG_DIST) & (Y < Y[i] +
                                          INIT_ANG_DIST) & (MAG < MAG[i] - RANGE_MAG)
        X_, Y_, MAG_ = X[cond], Y[cond], MAG[cond]
        dist2 = np.sort((X[i] - X_)**2 + (Y[i] - Y_)**2)
        if len(dist2) >= 2:
            min_dist[i] = sqrt(dist2[1])
    return min_dist


def create_pst_based_on_multipar(phot_file, mag_low, band, pst_file, nmax_psf, INIT_ANG_DIST, MAX_MSKY):
    """This function creates a new paremeter based on sharpness, groundness
    and sroundness, and uses this parameter to select stars to PSF.

    Parameters
    ----------
    infile : str
    Input file from detections.
    mag_low : float
    Limiting magnitude from the brightest star to select fainter candidate
    to PSF.
    band : str
    Name of the band.
    pst_file : str
    File name of DAOphot's task pstselect.
    nmax_psf : int
    Maximum number of candidates to PSF.
    """

    tilename = pst_file[0:13]

    ID, XCENTER, YCENTER, SHARPNESS, SROUND, GROUND, MAG = np.loadtxt(
        phot_file, usecols=(0, 1, 2, 4, 5, 6, 10), unpack=True)

    mean_sharp = np.mean(SHARPNESS)
    median_sharp = np.median(SHARPNESS)
    mean_sround = np.mean(SROUND)
    mean_ground = np.mean(GROUND)
    std_sharp = np.std(SHARPNESS)
    std_sround = np.std(SROUND)
    std_ground = np.std(GROUND)
    sharp_norm = (SHARPNESS - median_sharp)/std_sharp
    sround_norm = (SROUND - mean_sround)/std_sround
    ground_norm = (GROUND - mean_ground)/std_ground

    newpar = np.sqrt((sharp_norm**2 + sround_norm**2 + ground_norm**2) / 3)

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(
        nrows=1, ncols=4, sharey=False, figsize=(18, 6))
    ax0.hist(SHARPNESS, bins=30, range=(np.min(SHARPNESS), np.max(
        SHARPNESS)), histtype='step', label='Median={:.3f}'.format(median_sharp))
    ax0.set_xlabel('Sharpness')
    ax0.set_ylabel('counts')
    ax0.grid()
    ax0.legend()
    ax1.hist(SROUND, bins=30, range=(np.min(SROUND), np.max(SROUND)),
             histtype='step', label='Mean={:.3f}'.format(mean_sround))
    ax1.set_xlabel('Sround')
    ax1.grid()
    ax1.legend()
    ax2.hist(GROUND, bins=30, range=(np.min(GROUND), np.max(GROUND)),
             histtype='step', label='Mean={:.3f}'.format(mean_ground))
    ax2.set_xlabel('ground')
    ax2.grid()
    ax2.legend()
    ax3.hist(MAG, bins=30, range=(np.min(MAG), np.max(MAG)),
             histtype='step', label='Mag')
    ax3.set_xlabel('Mag')
    ax3.legend()
    ax3.set_yscale('log')
    plt.savefig(image_name + '_hist_pars.png')
    plt.close()

    plt.hist(sharp_norm, bins=30, range=(-4, 4),
             histtype='step', label='sharp_norm', color='b')
    plt.hist(sround_norm, bins=30, range=(-4, 4),
             histtype='step', label='sround_norm', color='g')
    plt.hist(ground_norm, bins=30, range=(-4, 4),
             histtype='step', label='ground_norm', color='r')
    plt.hist(newpar, bins=30, range=(-4, 4), histtype='step',
             label='new_par', color='maroon')
    plt.legend()
    plt.savefig(image_name + '_newpar.png')
    plt.close()

    min_d = min_dist(XCENTER, YCENTER, MAG, 7., INIT_ANG_DIST)

    bright_star = (MAG < np.min(MAG) + mag_low) & (min_d > INIT_ANG_DIST)
    XCENTER, YCENTER, ID, SHARPNESS, SROUND, GROUND, MAG, sharp_norm, sround_norm, ground_norm, newpar = XCENTER[bright_star], YCENTER[bright_star], ID[bright_star], SHARPNESS[
        bright_star], SROUND[bright_star], GROUND[bright_star], MAG[bright_star], sharp_norm[bright_star], sround_norm[bright_star], ground_norm[bright_star], newpar[bright_star]

    idx_par_sort = np.argsort(newpar)

    f = open(image_name + '_newpar_band_' + band + '.dat', 'w')
    for i in idx_par_sort:
        print(int(ID[i]), XCENTER[i], YCENTER[i],
              sround_norm[i], ground_norm[i], newpar[i], file=f)
    f.close()

    with open(pst_file) as ff:
        line = ff.read().splitlines()
    ff.close()

    g = open(pst_file[0:29] + '_imp_pst', 'w')
    for i in range(65):
        print(line[i], file=g)

    idx_pst, xcenter, ycenter, mag, msky = np.loadtxt(pst_file, unpack=True)

    idx_pst = [int(i) for i in idx_pst]
    count = 0
    for i in range(len(idx_par_sort)):
        if count < nmax_psf:
            if idx_par_sort[i] in idx_pst:
                idx_ = np.where(idx_pst == idx_par_sort[i])[0][0]
                if msky[idx_] < MAX_MSKY:
                    print("{:<9d}{:<10.3f}{:<10.3f}{:<12.3f}{:<15.7f}".format(
                        idx_pst[idx_], xcenter[idx_], ycenter[idx_], mag[idx_], msky[idx_]), file=g)
                    count += 1
        else:
            break
    g.close()

    return pst_file[0:29] + '_imp_pst'


def PSF_phot(fits_image, psf_file):
    """ Runs PSF photometry on image.

    Parameters
    ----------
    fits_image : str
        File name of image.
    imp_pst_file : str
        File name of improved PSF candidates.
    """
    print('Start to run allstar on {}'.format(fits_image))
    iraf.allstar(image=fits_image+"[0]", photfile='default', psfimage=psf_file, cache='no', allstarfile='default',
                 rejfile='default', subimage='default', verify='no', verbose='no')
    print('Finished to run allstar on {}'.format(fits_image))
    iraf.pdump(infile=fits_image+'0.als.1', fields="ID,XCENTER,YCENTER,MAG,MERR,SHARPNESS,CHI",
               expr="(MAG != INDEF)&&(MERR != INDEF)", Stdout=fits_image+'_pre_cal.dat')
    print('PSF_phot run on {}'.format(fits_image))


image_name = glob.glob('../1_DES_images_and_cats/images/*.fits')
tile_name_band = [(i.split('/')[3]).split('.')[0] for i in image_name]
a = [os.popen('wc -l ../3_ZP_fitting/' + i + '._calibrated_with_wcs.dat').read() for i in tile_name_band]
n_total_stars = [int(i.split(' ')[0]) for i in a]
n_total_ins_stars = [int(0.15 * i) for i in n_total_stars]
psf_files = [i + '.fits0.psf.1.fits' for i in tile_name_band]

for ii, jj in enumerate(tile_name_band):
    for kk in range(10):
        ins_image = tile_name_band[ii] + '.fits0.add.' + str(kk+1) + '.fits'
        print('Inserting stars in {} for the {:d} time'.format(ins_image, kk))

        add_star(image_name[ii], "", psf_files[ii], ins_image, 19., 28, n_total_ins_stars[ii])

        # ZP = zp_from_file(tile_name_band[ii])

        detection(ins_image)

        det_file, coo_file = ins_image + '_det.dat', ins_image + '0.coo.1'

        phot_pdump_file = ins_image + '0.mag'

        ap_phot(ins_image, det_file)

        PSF_phot(ins_image, psf_files[ii])

subprocess.call(['speech-dispatcher'])
subprocess.call(['spd-say', '"Your process has finished"'])
