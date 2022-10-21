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

mag_lim = 25

tiles = ['DES0052-2623_grizY__calibrated_with_wcs.dat', 'DES0312-5457_grizY__calibrated_with_wcs.dat', 'DES0057-3332_grizY__calibrated_with_wcs.dat', 'DES0356-4957_grizY__calibrated_with_wcs.dat', 'DES0057-3415_grizY__calibrated_with_wcs.dat', 'DES0424-2124_grizY__calibrated_with_wcs.dat', 'DES0100-3332_grizY__calibrated_with_wcs.dat', 'DES0511-3957_grizY__calibrated_with_wcs.dat', 'DES0100-3415_grizY__calibrated_with_wcs.dat', 'DES0515-3957_grizY__calibrated_with_wcs.dat', 'DES0203-0333_grizY__calibrated_with_wcs.dat', 'DES0523-2415_grizY__calibrated_with_wcs.dat', 'DES0221-0958_grizY__calibrated_with_wcs.dat', 'DES2134-0041_grizY__calibrated_with_wcs.dat', 'DES0224-0958_grizY__calibrated_with_wcs.dat']

for i in tiles:
    RA, DEC, MAGG = np.loadtxt('join_cats/calibrated_tiles/' + i, usecols=(1,2,3), unpack=True)
    RA, DEC, MAGG = RA[np.abs(MAGG) < mag_lim], DEC[np.abs(MAGG) < mag_lim], MAGG[np.abs(MAGG) < mag_lim]

    hdu = fits.open('cats/'+i[:12]+'_gold_y6.fits', memmap=True)
    RA_DES = hdu[1].data.field('ra')
    DEC_DES = hdu[1].data.field('dec')
    MAGG_DES = hdu[1].data.field('wavg_mag_psf_g')
    hdu.close()

    RA_DES, DEC_DES, MAGG_DES = RA_DES[MAGG_DES < mag_lim], DEC_DES[MAGG_DES < mag_lim], MAGG_DES[MAGG_DES < mag_lim]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, tight_layout=True, figsize=(21, 7))

    H1, xedges, yedges = np.histogram2d(RA, DEC, bins=[50,50], range=[[min(RA), max(RA)],[min(DEC), max(DEC)]])
    H1 = np.fliplr(np.flipud(H1.T))
    H2, xedges, yedges = np.histogram2d(RA_DES, DEC_DES, bins=[50,50], range=[[min(RA), max(RA)],[min(DEC), max(DEC)]])
    H2 = np.fliplr(np.flipud(H2.T))

    vmax = max(np.max(H1), np.max(H2))
    ax1.imshow(H1, aspect='auto', interpolation='None', cmap=cmap, extent=[xedges[-1], xedges[0], yedges[0], yedges[-1]], vmin=0, vmax=vmax)
    ax1.set_xlabel('RA')
    ax1.set_ylabel('DEC')
    ax1.set_title('DAOPHOT')

    ax2.imshow(H2, aspect='auto', interpolation='None', cmap=cmap, extent=[xedges[-1], xedges[0], yedges[0], yedges[-1]], vmin=0, vmax=vmax)
    ax2.set_xlabel('RA')

    ax2.set_title('DES')

    H3 = H1 / H2

    h3 = ax3.imshow(H3, aspect='auto', interpolation='None', cmap=cmap, extent=[xedges[-1], xedges[0], yedges[0], yedges[-1]])
    ax3.set_xlabel('RA')
    ax3.set_title('DAO/DES')
    plt.suptitle(i[:14])
    fig.colorbar(h3,fraction=0.046, pad=0.01)
    plt.savefig(i[:14] + '_ra_dec_dao_DES.png')
    plt.close()
