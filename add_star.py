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

def sky_sigma(im_data):
    sky = np.median(im_data)
    sigma = np.std(im_data[im_data < sky-np.min(im_data)])
    return sky, sigma
    
def add_star(image, photfile, psfimage, addimage, minmag, maxmag, nstar):
    iraf.daopars.unlearn
    iraf.photpars.unlearn
    iraf.datapars.unlearn
    iraf.fitskypars.unlearn
    iraf.centerpars.unlearn
    iraf.findpars.unlearn
    iraf.photpars.zmag = 25.0
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
    iraf.datapars.sigma = 4.0  # de acordo com a imagem
    iraf.datapars.readnoi = 0.0
    iraf.daopars.psfrad = 20.0  # de acordo com imagem # raio da maior imagem
    iraf.datapars.epadu = 1.  # de acordo com imagem
    iraf.daopars.fitrad = 6.  # de acordo com a imagem #(~1.5 ou 2.5 fwhm)
    iraf.daopars.recenter = 'yes'
    iraf.daopars.fitsky = 'yes'
    iraf.daopars.groupsk = 'yes'
    iraf.daopars.sannulu = 60.0
    iraf.fitskypars.annulus = 60
    iraf.fitskypars.dannulu = 10
    iraf.fitskypars.skyvalu = 0.1
    iraf.photpars.apertur = 12.0
    iraf.centerpars.cbox = 3
    iraf.centerpars.cthresh = 3.0
    iraf.addstar(image=image, photfile=photfile, psfimage=psfimage, addimage="default", minmag=minmag, maxmag=maxmag, nstar= nstar, verify='no',verbose='no')
    # Imprime a confirmacao da conclusao da tarefa
    print('task addstar concluida em {}'.format(image))

files = np.unique(glob.glob('*.fits'))
images = glob.glob('*.fits')[0]
tiles = np.unique([i[0:13] for i in files])


# Maybe create a folder to each tile and write outcomes on them.
bands = ['g', 'r', 'i']
#tirar ponto zero da referencia list_gc_tiles
#zp_g =
#zp_r =

'''
for ii in files:
    for jj in ['g', 'r']:
        tilename = ii[0:13] + jj
        
        image_name = glob.glob(ii[0:13] + jj + '.fits')[0]
'''  
     
for _ in range(10):
    #add_star(image_name[0], '', image_name[0] + "0.psf.1.fits", tilename + '_add_stars.fits', 19.0, 26.0, 7000) 
    add_star('DES0052-2623_i.fits[0]', '', 'DES0052-2623_i.fits0.psf.1.fits', 'DES0052-2623_i_add_stars.fits', 19, 26, 7000)
         
import subprocess
subprocess.call(['speech-dispatcher'])
subprocess.call(['spd-say', '"your process has finished"']) 

#fazer com q esse codigo crie automaticamente 10 imagens com estrelas artificiais inseridas com certo intervalo de magnitude 


