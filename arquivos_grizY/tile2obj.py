import numpy as np
import astropy.io.fits as fits

obj, tile = np.loadtxt('list_gc_tiles', usecols=(0,1), dtype=str, unpack=True)

objs_un, idx = np.unique(obj, return_inverse=True)

for i in range(len(objs_un)):

    name_file = objs_un[i]
    
    id_g, ra_g, dec_g, mag_g, magerr_g, sharp_g, chi_g, id_r, ra_r, dec_r, mag_r, magerr_r, sharp_r, chi_r, id_i, ra_i, dec_i, mag_i, magerr_i, sharp_i, chi_i, id_z, ra_z, dec_z, mag_z, magerr_z, sharp_z, chi_z, id_y, ra_y, dec_y, mag_y, magerr_y, sharp_y, chi_y, tile_a = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

    print(name_file)

    for j in tile[idx == i]:
        print(j)
        id_g_, ra_g_, dec_g_, mag_g_, magerr_g_, sharp_g_, chi_g_, id_r_, ra_r_, dec_r_, mag_r_, magerr_r_, sharp_r_, chi_r_, id_i_, ra_i_, dec_i_, mag_i_, magerr_i_, sharp_i_, chi_i_, id_z_, ra_z_, dec_z_, mag_z_, magerr_z_, sharp_z_, chi_z_, id_y_, ra_y_, dec_y_, mag_y_, magerr_y_, sharp_y_, chi_y_ = np.loadtxt(j + '_grizY__calibrated_with_wcs.dat', unpack=True)

        id_g.extend(id_g_)
        ra_g.extend(ra_g_)
        dec_g.extend(dec_g_)
        mag_g.extend(mag_g_)
        magerr_g.extend(magerr_g_)
        sharp_g.extend(sharp_g_)
        chi_g.extend(chi_g_)
        id_r.extend(id_r_)
        ra_r.extend(ra_r_)
        dec_r.extend(dec_r_)
        mag_r.extend(mag_r_)
        magerr_r.extend(magerr_r_)
        sharp_r.extend(sharp_r_)
        chi_r.extend(chi_r_)
        id_i.extend(id_i_)
        ra_i.extend(ra_i_)
        dec_i.extend(dec_i_)
        mag_i.extend(mag_i_)
        magerr_i.extend(magerr_i_)
        sharp_i.extend(sharp_i_)
        chi_i.extend(chi_i_)
        id_z.extend(id_z_)
        ra_z.extend(ra_z_)
        dec_z.extend(dec_z_)
        mag_z.extend(mag_z_)
        magerr_z.extend(magerr_z_)
        sharp_z.extend(sharp_z_)
        chi_z.extend(chi_z_)
        id_y.extend(id_y_)
        ra_y.extend(ra_y_)
        dec_y.extend(dec_y_)
        mag_y.extend(mag_y_)
        magerr_y.extend(magerr_y_)
        sharp_y.extend(sharp_y_)
        chi_y.extend(chi_y_)
        tile_a.extend(np.repeat(j, len(id_g_)))

    col1 = fits.Column(name='id_g', format='I', array=id_g)
    col2 = fits.Column(name='ra_g', format='D', array=ra_g)
    col3 = fits.Column(name='dec_g', format='D', array=dec_g)
    col4 = fits.Column(name='mag_g', format='E', array=mag_g)
    col5 = fits.Column(name='magerr_g', format='E', array=magerr_g)
    col6 = fits.Column(name='sharp_g', format='E', array=sharp_g)
    col7 = fits.Column(name='chi_g', format='E', array=chi_g)

    col8 = fits.Column(name='id_r', format='I', array=id_r)
    col9 = fits.Column(name='ra_r', format='D', array=ra_r)
    col10 = fits.Column(name='dec_r', format='D', array=dec_r)
    col11 = fits.Column(name='mag_r', format='E', array=mag_r)
    col12 = fits.Column(name='magerr_r', format='E', array=magerr_r)
    col13 = fits.Column(name='sharp_r', format='E', array=sharp_r)
    col14 = fits.Column(name='chi_r', format='E', array=chi_r)

    col15 = fits.Column(name='id_i', format='I', array=id_i)
    col16 = fits.Column(name='ra_i', format='D', array=ra_i)
    col17 = fits.Column(name='dec_i', format='D', array=dec_i)
    col18 = fits.Column(name='mag_i', format='E', array=mag_i)
    col19 = fits.Column(name='magerr_i', format='E', array=magerr_i)
    col20 = fits.Column(name='sharp_i', format='E', array=sharp_i)
    col21 = fits.Column(name='chi_i', format='E', array=chi_i)

    col22 = fits.Column(name='id_z', format='I', array=id_z)
    col23 = fits.Column(name='ra_z', format='D', array=ra_z)
    col24 = fits.Column(name='dec_z', format='D', array=dec_z)
    col25 = fits.Column(name='mag_z', format='E', array=mag_z)
    col26 = fits.Column(name='magerr_z', format='E', array=magerr_z)
    col27 = fits.Column(name='sharp_z', format='E', array=sharp_z)
    col28 = fits.Column(name='chi_z', format='E', array=chi_z)

    col29 = fits.Column(name='id_y', format='I', array=id_y)
    col30 = fits.Column(name='ra_y', format='D', array=ra_y)
    col31 = fits.Column(name='dec_y', format='D', array=dec_y)
    col32 = fits.Column(name='mag_y', format='E', array=mag_y)
    col33 = fits.Column(name='magerr_y', format='E', array=magerr_y)
    col34 = fits.Column(name='sharp_y', format='E', array=sharp_y)
    col35 = fits.Column(name='chi_y', format='E', array=chi_y)
    col36 = fits.Column(name='tile', format='A12', array=tile_a)

    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33, col34, col35, col36])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(name_file + '.fits')
    print('File for {} object done.'.format(name_file))
