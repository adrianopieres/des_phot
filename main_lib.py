from __future__ import print_function
import numpy as np
from math import sqrt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
# from stsci.tools import capable
#capable.OF_GRAPHICS = True
import pyraf
from pyraf import iraf
from pyraf.iraf import digiphot, daophot, imarith
import astropy.io.fits as fits
from pyraf import gwm
import glob, os
import matplotlib.pyplot as plt

# conda create -n iraf_py36 python=3.6
# conda activate iraf_py36
# sudo apt install python3-pyraf
# conda install matplotlib
# iraf_env


def ap_phot(im_name, coord_file, max_n_psf):
    """_summary_

    Parameters
    ----------
    im_name : _type_
        _description_
    coord_file : _type_
        _description_
    det_file : _type_
        _description_
    max_n_psf : _type_
        _description_
    """
    iraf.phot(image=im_name+'[0]', coord=coord_file,
              output='default', verbose='no', verify='no')
    iraf.pdump(infile=im_name + '0.mag.1', fields="ID,XCENTER,YCENTER,MAG,MERR,SHARPNESS,CHI",
               expr="(MAG != INDEF)&&(MERR != INDEF)", Stdout=im_name + '0.mag')
    iraf.pstselect(image=im_name+'[0]', photfile='default',
                   pstfile='default', maxnpsf=max_n_psf, verbose='no', verify='no')
    print('Ap_phot run for {}'.format(im_name))


def detection(im_name):
    """Run script in the detections image

    Parameters
    ----------
    im_name : _type_
        _description_
    """
    iraf.daopars.unlearn
    iraf.photpars.unlearn
    iraf.datapars.unlearn
    iraf.fitskypars.unlearn
    iraf.centerpars.unlearn
    iraf.findpars.unlearn
    iraf.photpars.zmag = 30.0
    # Procurar no ds9 #Maximum equivalent exposure time (s)
    iraf.datapars.exposure = "EXPTIME"
    iraf.daopars.functio = "auto"
    iraf.datapars.gain = "GAIN"  # Maximum equivalent gain (e-/ADU)
    iraf.datapars.scale = 1
    iraf.datapars.datamin = -100.0
    iraf.datapars.datamax = 1500.0  # ou utilizar #1500
    iraf.daopars.saturated = 'yes'
    # de acordo com a imagem #desvio padra da #imexamin(comando 'r') no ds9
    iraf.datapars.fwhmpsf = 11.84
    iraf.datapars.sigma = 0.0275  # de acordo com a imagem
    iraf.datapars.readnoi = 0.0
    iraf.daopars.psfrad = 161.5  # de acordo com imagem # raio da maior imagem
    iraf.datapars.epadu = 1  # de acordo com imagem
    iraf.daopars.fitrad = 4  # de acordo com a imagem #(~1.5 ou 2.5 fwhm)
    iraf.daopars.recenter = 'yes'
    iraf.daopars.fitsky = 'yes'
    iraf.daopars.groupsk = 'yes'
    iraf.daopars.sannulu = 20.0
    iraf.fitskypars.annulus = 10
    iraf.fitskypars.dannulu = 6
    iraf.fitskypars.skyvalu = 0.004
    iraf.photpars.apertur = 6.0
    iraf.centerpars.cbox = 3
    iraf.centerpars.cthresh = 3.0

    iraf.daofind(image=im_name + '[0]', threshold=3.00,
                 output='default', verbose='no', verify='no')
    iraf.pdump(infile=im_name + '0.coo.1', fields="ID, XCENTER, YCENTER, MAG, SHARPNESS, SROUND, GROUND",
               expr="(MAG != INDEF)&&(SHARPNESS != INDEF)&&(SROUND != INDEF)&&(GROUND != INDEF)", Stdout=im_name + '0_det.dat')


def matching_stars(infile_phot, infile_DES):
    """_summary_

    Parameters
    ----------
    infile_phot : _type_
        _description_
    infile_DES : _type_
        _description_
    """
    RA, DEC, MAG_I, MAGI_ERR = np.loadtxt(
        infile_phot, usecols=(1, 2, 3, 4), unpack=True)
    RA_DES, DEC_DES, MAGG_DES, MAGG_ERR_DES, MAGR_DES, MAGR_ERR_DES, MAGI_DES, MAGI_ERR_DES, MAGZ_DES, MAGZ_DES_ERR, MAGY_DES, MAGY_ERR_DES = np.loadtxt(
        infile_DES, usecols=(2, 3, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17), delimiter=',', skiprows=1, unpack=True)

    # Agora configurando cada catalogo
    C_DAO = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)
    C_DES = SkyCoord(ra=RA_DES*u.degree, dec=DEC_DES*u.degree)

    idx_dao, idx_des, d2d, d3d = C_DES.search_around_sky(C_DAO, 1*u.arcsec)

    # agora escrevendo ambas as coordenadas em ambos os arquivos:
    for i, j in zip(idx_dao, idx_des):
        print(C_DAO[i].ra.deg, C_DAO[i].dec.deg, MAG_I[i], MAGI_ERR[i],
              C_DES[j].ra.deg, C_DES[j].dec.deg, MAGI_DES[j], MAGI_ERR_DES[j], file=f)
    f.close()


def calc_ZP(infile, infile_DES, band, mag_lim_sat_DES):
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

    # read files...
    # excluir estrelas saturadas com mag < mag_lim_sat:
    idx_des_faint = [i for i in idx_des if MAG_DES[j] > 17.]
    idx_dao_faint = [(i, j) for i, j in zip(idx_des) if MAG_DES[j] > 17.]

    # plot Mag fotometria x Mag DES:

    best_fit = np.polyfit(
        MAG[idx_dao_faint], MAG_DES[idx_des_faint], deg=1, w=1/MAG_ERR[idx_dao_faint]**2)

    plt.scatter(MAG[idx_dao_faint], MAG_DES[idx_des_faint])
    plt.xlabel('MAG in i band')
    plt.legend()
    plt.xlabel('Magnitude in daophot')
    plt.ylabel('Magnitude DES')
    plt.title('Best-Fitting for {} band: {:.4f}, {:.4f}'.format(band,
              best_fit[0], best_fit[1]))
    plt.savefig(infile[0:13] + '_ZP_' + band + '.png')
    plt.close()
    return best_fit[0]


def min_dist(X, Y, MAG, RANGE_MAG, INIT_ANG_DIST=20.):
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


def create_pst_based_on_multipar(phot_file, mag_low, band, pst_file, nmax_psf, INIT_ANG_DIST):
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
    print(pst_file)

    ID, XCENTER, YCENTER, SHARPNESS, SROUND, GROUND, MAG = np.loadtxt(
        phot_file, usecols=(0, 1, 2, 4, 5, 6, 9), unpack=True)

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

    # Are we selecting good stars this way?
    # Maybe create a parameter with restricted stars from DES?
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
    plt.savefig(tilename + '_hist_pars.png')
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
    plt.savefig(tilename + '_newpar.png')
    plt.close()

    min_dist(XCENTER, YCENTER, MAG, 5., INIT_ANG_DIST)

    bright_star = (MAG < np.min(MAG) + mag_low) & (min_dist < INIT_ANG_DIST)
    XCENTER, YCENTER, ID, SHARPNESS, SROUND, GROUND, MAG, sharp_norm, sround_norm, ground_norm, newpar = XCENTER[bright_star], YCENTER[bright_star], ID[bright_star], SHARPNESS[
        bright_star], SROUND[bright_star], GROUND[bright_star], MAG[bright_star], sharp_norm[bright_star], sround_norm[bright_star], ground_norm[bright_star], newpar[bright_star]

    idx_par_sort = np.argsort(newpar)
    id_sort = [ID[i] for i in idx_par_sort]

    f = open(tilename + '_newpar_band_' + band + '.dat', 'w')
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
    print('PST_FILE == ', pst_file)
    idx_pst, xcenter, ycenter, mag, msky = np.loadtxt(
        pst_file, unpack=True)
    idx_pst = [int(i) for i in idx_pst]

    count = 0
    while count < nmax_psf:
        for i in range(len(idx_par_sort)):
            if idx_par_sort[i] in idx_pst:
                idx_ = np.where(idx_pst == idx_par_sort[i])[0][0]
    # inseri < para indexar a esquerda e parece funcionar
                print("{:<9d}{:<10.3f}{:<10.3f}{:<12.3f}{:<15.7f}".format(
                    idx_pst[idx_], xcenter[idx_], ycenter[idx_], mag[idx_], msky[idx_]), file=g)
                count += 1
    g.close()
    return pst_file[0:29] + '_imp_pst'


def PSF_phot(fits_image, imp_pst_file):
    # Run on iraf:
    # psf_filename = fits_image[0:13] + '_psf.fits'
    iraf.psf(image=fits_image + "[0]", photfile='default', pstfile=imp_pst_file,
             psfimage='default', opstfile='default', groupfile='default', interactive='no',
             verify='no', verbose='no')
    iraf.seepsf(psfimage=fits_image + "0.fits.psf",
                image=fits_image + '_seepsf.fits')
    iraf.allstar(image=fits_image+"[0]", photfile='default', psfimage='default',
                 cache='no', allstarfile='default', rejfile='default',
                 subimage='default', verify='no', verbose='no')
    iraf.pdump(infile=image+'0.als.1', fields="ID,XCENTER,YCENTER,MAG,MERR,SHARPNESS,CHI",
               expr="(MAG != INDEF)&&(MERR != INDEF)", Stdout=image+'_pre_cal.dat')

    print('PSF_phot run on {}'.format(fits_image))


def wcs(im_name, infile, outfile):
    iraf.wcsctran(input=infile, output=outfile, image=im_name,
                  inwcs='logical', outwcs='world', col="2 3")
    print('wcs run on {}'.format(outfile))


files = glob.glob('*.fits')
det_images = glob.glob('*_det.fits')
tiles = [i[0:13] for i in det_images]
# Maybe create a folder to each tile and write outcomes on them.
bands = ['g', 'r', 'i', 'z', 'Y']

# According to DES DR2 paper this are the limits for saturation:
# lim_mag = [15.2, 15.7, 15.8, 15.5, 13.6]
# We are setting 2 magnitudes fainter than saturation limit:
lim_mag = [17.2, 17.7, 17.8, 17.5, 15.6]

for ii in det_images:
    detection(ii)
    det_file, coo_file = ii + '0_det.dat', ii + '0.coo.1'
    for jj in bands:
        tilename = ii[0:13] + jj
        image_name = glob.glob(ii[0:13] + jj + '.fits')[0]
        phot_pdump_file = image_name + '0.mag'
        
        ap_phot(image_name, coo_file, 200)
        
        pst_file = glob.glob(ii[0:13] + jj + '*.pst.1')[0]
        
        os.system('join --nocheck-order ' + det_file  + ' ' + phot_pdump_file + ' > ' + tilename + '_parsfile.dat')

        imp_pst_file = create_pst_based_on_multipar(tilename + '_parsfile.dat',
                                                    3., jj, pst_file, 50, 50)
        all_star_file_flat = PSF_phot(image_name, imp_pst_file)

        wcs(image_name, all_star_file_flat, tilename + '_wcs_not_cal.dat')

        ZP = calc_ZP(tilename + '_wcs_not_cal.dat', tilename +
                     '_' + jj + '.csv', jj, mag_lim_sat_DES)

import subprocess
subprocess.call(['speech-dispatcher'])
subprocess.call(['spd-say', '"your process has finished"'])
