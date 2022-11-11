#hist2d com RA DEC, bins = 100
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
import astropy.io.fits as fits
import matplotlib as mpl
mpl.rcParams['legend.numpoints'] = 1
cmap = mpl.cm.get_cmap("inferno_r")
cmap2 = mpl.cm.get_cmap("inferno")
from matplotlib.colors import LogNorm

mag_lim = 25

tiles = ['DES0052-2623_grizY__calibrated_with_wcs.dat', 'DES0312-5457_grizY__calibrated_with_wcs.dat', 'DES0057-3332_grizY__calibrated_with_wcs.dat', 'DES0356-4957_grizY__calibrated_with_wcs.dat', 'DES0057-3415_grizY__calibrated_with_wcs.dat', 'DES0424-2124_grizY__calibrated_with_wcs.dat', 'DES0100-3332_grizY__calibrated_with_wcs.dat', 'DES0511-3957_grizY__calibrated_with_wcs.dat', 'DES0100-3415_grizY__calibrated_with_wcs.dat', 'DES0515-3957_grizY__calibrated_with_wcs.dat', 'DES0203-0333_grizY__calibrated_with_wcs.dat', 'DES0523-2415_grizY__calibrated_with_wcs.dat', 'DES0221-0958_grizY__calibrated_with_wcs.dat', 'DES2134-0041_grizY__calibrated_with_wcs.dat', 'DES0224-0958_grizY__calibrated_with_wcs.dat']

for i in tiles:
    RA, DEC, MAGG, MAGR, MAGI, MAGZ, MAGY = np.loadtxt('join_cats/calibrated_tiles/' + i, usecols=(1,2,3,10,17,24,31), unpack=True)
    RA, DEC, MAGG, MAGR, MAGI, MAGZ, MAGY = RA[np.abs(MAGG) < mag_lim], DEC[np.abs(MAGG) < mag_lim], MAGG[np.abs(MAGG) < mag_lim], MAGR[np.abs(MAGG) < mag_lim], MAGI[np.abs(MAGG) < mag_lim], MAGZ[np.abs(MAGG) < mag_lim], MAGY[np.abs(MAGG) < mag_lim]
    GR = MAGG-MAGR
    RI = MAGR-MAGI
    IZ = MAGI-MAGZ
    ZY = MAGZ-MAGY

    hdu = fits.open('cats/'+i[:12]+'_gold_y6.fits', memmap=True)
    RA_DES = hdu[1].data.field('ra')
    DEC_DES = hdu[1].data.field('dec')
    MAGG_DES = hdu[1].data.field('wavg_mag_psf_g')
    MAGR_DES = hdu[1].data.field('wavg_mag_psf_r')
    MAGI_DES = hdu[1].data.field('wavg_mag_psf_i')
    MAGZ_DES = hdu[1].data.field('wavg_mag_psf_z')
    MAGY_DES = hdu[1].data.field('wavg_mag_psf_Y')
    hdu.close()

    RA_DES, DEC_DES, MAGG_DES, MAGR_DES, MAGI_DES, MAGZ_DES, MAGY_DES = RA_DES[MAGG_DES < mag_lim], DEC_DES[MAGG_DES < mag_lim], MAGG_DES[MAGG_DES < mag_lim], MAGR_DES[MAGG_DES < mag_lim], MAGI_DES[MAGG_DES < mag_lim], MAGZ_DES[MAGG_DES < mag_lim], MAGY_DES[MAGG_DES < mag_lim]
    GR_DES = MAGG_DES - MAGR_DES
    RI_DES = MAGR_DES - MAGI_DES
    IZ_DES = MAGI_DES - MAGZ_DES
    ZY_DES = MAGZ_DES - MAGY_DES

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, tight_layout=False, figsize=(21, 7))

    H1, xedges, yedges = np.histogram2d(RA, DEC, bins=[50,50], range=[[min(RA), max(RA)],[min(DEC), max(DEC)]])
    H1 = np.fliplr(np.flipud(H1.T))
    H2, xedges, yedges = np.histogram2d(RA_DES, DEC_DES, bins=[50,50], range=[[min(RA), max(RA)],[min(DEC), max(DEC)]])
    H2 = np.fliplr(np.flipud(H2.T))

    vmax = max(np.max(H1), np.max(H2))
    h1 = ax1.imshow(H1, aspect='auto', interpolation='None', cmap=cmap, extent=[xedges[-1], xedges[0], yedges[0], yedges[-1]], vmin=0, vmax=vmax)
    ax1.set_xlabel('RA')
    ax1.set_ylabel('DEC')
    ax1.set_title('DAOPHOT')
    # plt.colorbar(h1, fraction=0.046, pad=0.5)
    cbaxes = fig.add_axes([0.355, 0.11, 0.015, 0.77])
    cbar = fig.colorbar(h1, cax=cbaxes, cmap=cmap, orientation='vertical')

    h2 = ax2.imshow(H2, aspect='auto', interpolation='None', cmap=cmap, extent=[xedges[-1], xedges[0], yedges[0], yedges[-1]], vmin=0, vmax=vmax)
    ax2.set_xlabel('RA')

    ax2.set_title('DES')
    cbaxes = fig.add_axes([0.632, 0.11, 0.015, 0.77])
    cbar = fig.colorbar(h2, cax=cbaxes, cmap=cmap, orientation='vertical')

    H3 = H1 / H2

    h3 = ax3.imshow(H3, aspect='auto', interpolation='None', cmap=cmap, extent=[xedges[-1], xedges[0], yedges[0], yedges[-1]])
    ax3.set_xlabel('RA')
    ax3.set_title('DAO/DES')
    plt.suptitle(i[:14])
    cbaxes = fig.add_axes([0.905, 0.11, 0.015, 0.77])
    cbar = fig.colorbar(h3, cax=cbaxes, cmap=cmap, orientation='vertical')
    plt.savefig(i[:14] + '_ra_dec_dao_DES.png')
    plt.close()
    
    #################################################
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, tight_layout=False, figsize=(15, 7))

    H1, xedges, yedges = np.histogram2d(GR, RI, bins=[50,50], range=[[-0.5, 2.0],[-0.5, 2.0]])
    H1 = np.flipud(H1.T)
    H1[H1 < 0.01] = 0.01
    H2, xedges, yedges = np.histogram2d(GR_DES, RI_DES, bins=[50,50], range=[[-0.5, 2.0],[-0.5, 2.0]])
    H2 = np.flipud(H2.T)
    H2[H2 < 0.01] = 0.01

    vmax = max(np.max(H1), np.max(H2))
    h1 = ax1.imshow(H1, aspect='auto', interpolation='None', cmap=cmap2, extent=[-0.5, 2.0, -0.5, 2.0], norm=LogNorm(vmin=0.01, vmax=vmax))
    ax1.set_xlabel('g-r')
    ax1.set_ylabel('r-i')
    ax1.set_title('DAOPHOT')

    h2 = ax2.imshow(H2, aspect='auto', interpolation='None', cmap=cmap2, extent=[-0.5, 2.0, -0.5, 2.0], norm=LogNorm(vmin=0.01, vmax=vmax))
    ax2.set_xlabel('g-r')
    ax2.set_ylabel('r-i')
    ax2.set_title('DES')
        
    cbaxes = fig.add_axes([0.905, 0.11, 0.015, 0.77])
    cbar = fig.colorbar(h2, cax=cbaxes, cmap=cmap2, orientation='vertical')

    plt.savefig(i[:14] + '_gri_color-color_dao_des.png')
    plt.close()
    #################################################
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, tight_layout=False, figsize=(15, 7))

    H1, xedges, yedges = np.histogram2d(IZ, ZY, bins=[50,50], range=[[-0.5, 2.0],[-0.5, 2.0]])
    H1 = np.flipud(H1.T)
    H1[H1 < 0.01] = 0.01
    H2, xedges, yedges = np.histogram2d(IZ_DES, ZY_DES, bins=[50,50], range=[[-0.5, 2.0],[-0.5, 2.0]])
    H2 = np.flipud(H2.T)
    H2[H2 < 0.01] = 0.01

    vmax = max(np.max(H1), np.max(H2))
    h1 = ax1.imshow(H1, aspect='auto', interpolation='None', cmap=cmap2, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], norm=LogNorm(vmin=0.01, vmax=vmax))
    ax1.set_xlabel('i-z')
    ax1.set_ylabel('z-Y')
    ax1.set_title('DAOPHOT')

    h2 = ax2.imshow(H2, aspect='auto', interpolation='None', cmap=cmap2, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], norm=LogNorm(vmin=0.01, vmax=vmax))
    ax2.set_xlabel('i-z')
    ax2.set_ylabel('z-Y')
    ax2.set_title('DES')
        
    cbaxes = fig.add_axes([0.905, 0.11, 0.015, 0.77])
    cbar = fig.colorbar(h2, cax=cbaxes, cmap=cmap2, orientation='vertical')

    plt.savefig(i[:14] + '_izY_color-color_dao_des.png')
    plt.close()
