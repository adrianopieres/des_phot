import numpy as np
import astropy.io.fits as fits

obj, tile = np.loadtxt('list_gc_tiles', usecols=(0, 1), dtype=str, unpack=True)

ramin, ramax, decmin, decmax = np.loadtxt('list_gc_tiles', usecols=(2, 3, 4, 5), unpack=True)

objs_un, idx = np.unique(obj, return_inverse=True)

for i in range(len(objs_un)):

    name_file = objs_un[i]

    id_g, x_g, y_g, ra_g, dec_g, mag_g, magerr_g, sharp_g, chi_g, id_r, x_r, y_r, ra_r, dec_r, mag_r, magerr_r, sharp_r, chi_r, id_i, x_i, y_i, ra_i, dec_i, mag_i, magerr_i, sharp_i, chi_i, id_z, x_z, y_z, ra_z, dec_z, mag_z, magerr_z, sharp_z, chi_z, id_y, x_y, y_y, ra_y, dec_y, mag_y, magerr_y, sharp_y, chi_y, tile_a = [
    ], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

    print(name_file)

    for j in tile[idx == i]:
        idx_ra = np.where(tile == j)[0]
        print(j, idx_ra)
        print('../4_join_bands/' + j + '_grizY_calibrated.dat')
        id_g_, x_g_, y_g_, ra_g_, dec_g_, mag_g_, magerr_g_, sharp_g_, chi_g_, id_r_, x_r_, y_r_, ra_r_, dec_r_, mag_r_, magerr_r_, sharp_r_, chi_r_, id_i_, x_i_, y_i_, ra_i_, dec_i_, mag_i_, magerr_i_, sharp_i_, chi_i_, id_z_, x_z_, y_z_, ra_z_, dec_z_, mag_z_, magerr_z_, sharp_z_, chi_z_, id_y_, x_y_, y_y_, ra_y_, dec_y_, mag_y_, magerr_y_, sharp_y_, chi_y_ = np.loadtxt('../4_join_bands/' + j + '_grizY_calibrated.dat', unpack=True)

        cond = ((ra_g_ > ramin[idx_ra])&(ra_g_ < ramax[idx_ra])&(dec_g_ > decmin[idx_ra])&(dec_g_ < decmax[idx_ra]))|((ra_r_ > ramin[idx_ra])&(ra_r_ < ramax[idx_ra])&(dec_r_ > decmin[idx_ra])&(dec_r_ < decmax[idx_ra]))|((ra_i_ > ramin[idx_ra])&(ra_i_ < ramax[idx_ra])&(dec_i_ > decmin[idx_ra])&(dec_i_ < decmax[idx_ra]))|((ra_z_ > ramin[idx_ra])&(ra_z_ < ramax[idx_ra])&(dec_z_ > decmin[idx_ra])&(dec_z_ < decmax[idx_ra]))|((ra_y_ > ramin[idx_ra])&(ra_y_ < ramax[idx_ra])&(dec_y_ > decmin[idx_ra])&(dec_y_ < decmax[idx_ra]))
        print(len(id_g_), len(id_g_[cond]))
        id_g.extend(id_g_[cond])
        x_g.extend(x_g_[cond])
        y_g.extend(y_g_[cond])
        ra_g.extend(ra_g_[cond])
        dec_g.extend(dec_g_[cond])
        mag_g.extend(mag_g_[cond])
        magerr_g.extend(magerr_g_[cond])
        sharp_g.extend(sharp_g_[cond])
        chi_g.extend(chi_g_[cond])
        id_r.extend(id_r_[cond])
        x_r.extend(x_r_[cond])
        y_r.extend(y_r_[cond])
        ra_r.extend(ra_r_[cond])
        dec_r.extend(dec_r_[cond])
        mag_r.extend(mag_r_[cond])
        magerr_r.extend(magerr_r_[cond])
        sharp_r.extend(sharp_r_[cond])
        chi_r.extend(chi_r_[cond])
        id_i.extend(id_i_[cond])
        x_i.extend(x_i_[cond])
        y_i.extend(y_i_[cond])
        ra_i.extend(ra_i_[cond])
        dec_i.extend(dec_i_[cond])
        mag_i.extend(mag_i_[cond])
        magerr_i.extend(magerr_i_[cond])
        sharp_i.extend(sharp_i_[cond])
        chi_i.extend(chi_i_[cond])
        id_z.extend(id_z_[cond])
        x_z.extend(x_z_[cond])
        y_z.extend(y_z_[cond])
        ra_z.extend(ra_z_[cond])
        dec_z.extend(dec_z_[cond])
        mag_z.extend(mag_z_[cond])
        magerr_z.extend(magerr_z_[cond])
        sharp_z.extend(sharp_z_[cond])
        chi_z.extend(chi_z_[cond])
        id_y.extend(id_y_[cond])
        x_y.extend(x_y_[cond])
        y_y.extend(y_y_[cond])
        ra_y.extend(ra_y_[cond])
        dec_y.extend(dec_y_[cond])
        mag_y.extend(mag_y_[cond])
        magerr_y.extend(magerr_y_[cond])
        sharp_y.extend(sharp_y_[cond])
        chi_y.extend(chi_y_[cond])
        tile_a.extend(np.repeat(j, len(id_g_[cond])))

    col1 = fits.Column(name='id_g', format='I', array=id_g)
    col2 = fits.Column(name='x_g', format='E', array=x_g)
    col3 = fits.Column(name='y_g', format='E', array=y_g)
    col4 = fits.Column(name='ra_g', format='D', array=ra_g)
    col5 = fits.Column(name='dec_g', format='D', array=dec_g)
    col6 = fits.Column(name='mag_g', format='E', array=mag_g)
    col7 = fits.Column(name='magerr_g', format='E', array=magerr_g)
    col8 = fits.Column(name='sharp_g', format='E', array=sharp_g)
    col9 = fits.Column(name='chi_g', format='E', array=chi_g)

    col10 = fits.Column(name='id_r', format='I', array=id_r)
    col11 = fits.Column(name='x_r', format='E', array=x_r)
    col12 = fits.Column(name='y_r', format='E', array=y_r)
    col13 = fits.Column(name='ra_r', format='D', array=ra_r)
    col14 = fits.Column(name='dec_r', format='D', array=dec_r)
    col15 = fits.Column(name='mag_r', format='E', array=mag_r)
    col16 = fits.Column(name='magerr_r', format='E', array=magerr_r)
    col17 = fits.Column(name='sharp_r', format='E', array=sharp_r)
    col18 = fits.Column(name='chi_r', format='E', array=chi_r)

    col19 = fits.Column(name='id_i', format='I', array=id_i)
    col20 = fits.Column(name='x_i', format='E', array=x_i)
    col21 = fits.Column(name='y_i', format='E', array=y_i)
    col22 = fits.Column(name='ra_i', format='D', array=ra_i)
    col23 = fits.Column(name='dec_i', format='D', array=dec_i)
    col24 = fits.Column(name='mag_i', format='E', array=mag_i)
    col25 = fits.Column(name='magerr_i', format='E', array=magerr_i)
    col26 = fits.Column(name='sharp_i', format='E', array=sharp_i)
    col27 = fits.Column(name='chi_i', format='E', array=chi_i)

    col28 = fits.Column(name='id_z', format='I', array=id_z)
    col29 = fits.Column(name='x_z', format='E', array=x_z)
    col30 = fits.Column(name='y_z', format='E', array=y_z)
    col31 = fits.Column(name='ra_z', format='D', array=ra_z)
    col32 = fits.Column(name='dec_z', format='D', array=dec_z)
    col33 = fits.Column(name='mag_z', format='E', array=mag_z)
    col34 = fits.Column(name='magerr_z', format='E', array=magerr_z)
    col35 = fits.Column(name='sharp_z', format='E', array=sharp_z)
    col36 = fits.Column(name='chi_z', format='E', array=chi_z)

    col37 = fits.Column(name='id_y', format='I', array=id_y)
    col38 = fits.Column(name='x_y', format='E', array=x_y)
    col39 = fits.Column(name='y_y', format='E', array=y_y)
    col40 = fits.Column(name='ra_y', format='D', array=ra_y)
    col41 = fits.Column(name='dec_y', format='D', array=dec_y)
    col42 = fits.Column(name='mag_y', format='E', array=mag_y)
    col43 = fits.Column(name='magerr_y', format='E', array=magerr_y)
    col44 = fits.Column(name='sharp_y', format='E', array=sharp_y)
    col45 = fits.Column(name='chi_y', format='E', array=chi_y)
    col46 = fits.Column(name='tile', format='A12', array=tile_a)

    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17,
                        col18, col19, col20, col21, col22, col23, col24, col25, col26, col27, col28, col29, col30, col31, col32, col33,
                        col34, col35, col36, col37, col38, col39, col40, col41, col42, col43, col44, col45, col46])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto('objects/' + name_file + '.fits', overwrite=True)
    print('File for {} object done.'.format(name_file))