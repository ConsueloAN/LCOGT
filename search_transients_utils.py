# Tiene todas las funciones utiles para search_transients.py

import numpy as np
import fitsio
import copy
from astropy.io import fits
from scipy import signal
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
from matplotlib import rcParams
from astropy.table import Table
from matplotlib.patches import Ellipse
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_skycoord
from astropy.wcs import WCS
from astropy import units as u
import math
import sep
import csv
import os
import numpy as np


def get_CC(image, skymapper, filtro, mag_max=16., plot=False):
	""" This function gets calibration constant to calculate flux by matching with skymapper catalog """

	# para calcular constante de calibracion CC:
	data = fitsio.read(image)
	bkg = sep.Background(data)
	data_sub = data - bkg

	objects = sep.extract(data_sub, 1.9, err=bkg.globalrms) # deteccion de objetos
	
	flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'], 8, err=bkg.globalrms, gain=1.0) # flujo de objetos

	# conversion a WSC 
	header = fits.getheader(image)
	w = WCS(header)
	wx, wy = w.wcs_pix2world(objects['x'], objects['y'], 1)

	cond = (flux >= 0) & (flux-fluxerr >= 0) # condiciones para log
	tab_SEP = Table()
	tab_SEP['inds_SEP'] = range(len(flux[cond]))
	tab_SEP['ra'] = wx[cond]
	tab_SEP['dec'] = wy[cond]
	tab_SEP['flux'] = flux[cond]

	flux = np.array(flux)
	mag_nc = []
	for i in tab_SEP['flux']:
		mag_nc += [-2.5 * math.log10(i)]
	tab_SEP['mag_nc'] = mag_nc

	skymapper = ascii.read(skymapper, format = 'tab') 
	skymapper_coords = SkyCoord(skymapper['raj2000'], skymapper['dej2000'], unit = 'deg')
	SEP_coords = SkyCoord(tab_SEP['ra'], tab_SEP['dec'], unit = 'deg')

	idx, sep2D, sep3D = SEP_coords.match_to_catalog_sky(skymapper_coords, nthneighbor = 1)

	tab_SEP['ind_skymapper'] = idx
	tab_SEP['sep2D'] = sep2D
	sep2D = np.array(sep2D)

	inds_zero_sep = np.where(sep2D <= 0.0005)[0]     # 2 arcsec
	tab_zero_sep = tab_SEP[inds_zero_sep]

	tab_zero_sep['mag_nc'] = tab_SEP['mag_nc'][tab_zero_sep['inds_SEP']]
	tab_zero_sep['mag_skymapper'] = skymapper[filtro][tab_zero_sep['ind_skymapper']]
	tab_zero_sep['dif_magskymapper_magnc'] = tab_zero_sep['mag_skymapper'] - tab_zero_sep['mag_nc'] 

	#condicion para eliminar peak en 0 & estrellas brillantes saturadas & alta dispersion en magnitudes grandes
	cond = (tab_zero_sep['dif_magskymapper_magnc'] >= 5.) & (tab_zero_sep['mag_skymapper'] <= mag_max)
	CC = np.median(tab_zero_sep['dif_magskymapper_magnc'][cond])

	# plot magnitud calibrada con CC vs magnitud skymapper
	mag_c = []
	for row in tab_SEP:
		flux_aux = row['flux']
		mag_c += [CC -2.5 * math.log10(flux_aux)]
	tab_SEP['mag_c'] = mag_c
	tab_zero_sep['mag_c'] = tab_SEP['mag_c'][tab_zero_sep['inds_SEP']]
	
	if plot:
		plt.plot(tab_zero_sep['mag_c'], tab_zero_sep['mag_skymapper'], 'om')
		plt.plot(range(10, 23), range(10, 23), 'y')
		plt.title("mag SEP vs mag skymapper")
		plt.ylabel('mag_skymapper')
		plt.xlabel('mag_SEP')
		plt.show()

		plt.plot(tab_zero_sep['mag_c'], tab_zero_sep['mag_c'] - tab_zero_sep['mag_skymapper'], 'om')
		plt.plot(range(10, 23), np.zeros(len(range(10, 23))), 'y')
		plt.title("mag SEP vs diff mag")
		plt.ylabel('mag_SEP - mag_skymapper')
		plt.xlabel('mag_SEP')
		plt.show()

	return CC  #constante de calibracion de magnitud

def get_flux(mag, CC):
	"""This function gets flux from a magnitude and CC"""

	# mag = -2.5*math.log10(flux) + CC
	flux = 10**(-(mag - CC)/2.5)
	return flux

def get_star_norm(size, FWHM):
	"""This function gets a normalized gaussian matrix"""

	std = FWHM/2.35482
	get_star_1d = signal.gaussian(size, std=std).reshape(size, 1)
	get_star_2d = np.outer(get_star_1d, get_star_1d)
	return get_star_2d

def get_star(star_norm, flux):
	"""This function gets a gaussian matrix with a determined flux"""

	# valor A que debo multiplicar a la matriz para que la suma de elementos sea igual al flujo:
	A = flux / np.sum(star_norm)
	star = star_norm * A
	return star

def get_mag_rec(image_diff, CC, circle_radio, x, y):
	"""This function gets the magnitude of a star at position x, y of a differentiated image"""

	data = fitsio.read(image_diff)
    # sustraccion background
	bkg = sep.Background(data)     # estima el background
	data_sub = data - bkg    # resta el background estimado
	objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)    # detecta objetos
	x = [x]
	y = [y]
	flux_s, fluxerr_s, flag_s = sep.sum_circle(data_sub, x, y, circle_radio, err=bkg.globalrms, gain=1.0)

	flux_s = np.where(flux_s>0, flux_s, 1)
	mas = flux_s + fluxerr_s
	menos = flux_s - fluxerr_s
	mas = np.where(mas>0, mas, 1)
	menos = np.where(menos>0, menos,1)
	mag = CC -2.5 * math.log10(flux_s)
	mag_m = CC -2.5 * math.log10(mas)    # porque mas flujo es menos magnitude
	mag_p = CC -2.5 * math.log10(menos)

	err_mag_m = mag - mag_m
	err_mag_p = mag_p - mag

	SN = float(flux_s/fluxerr_s)
	
	return (mag, err_mag_p, err_mag_m, SN)

def mag_rec(image, skymapper, out_image_star, filtro, mag, x, y, h_size, FWHM, circle_radio, image_template, image_diff):
	"""This function makes an artificial star of a given magnitude at the x, y position of an image, differentiates this 
	image (ToO) with a template from a different epoch and calculates the recovered magnitude of the artificial star in 
	the difference image ."""

	CC = get_CC(image, skymapper, filtro)
	flux = get_flux(mag, CC)
	star_norm = get_star_norm(2*h_size, FWHM)
	star = get_star(star_norm, flux)

	zeros = np.zeros((4096, 4096))
	y_menos = y - h_size
	y_mas = y + h_size
	x_menos = x - h_size
	x_mas = x + h_size
	zeros[y_menos:y_mas, x_menos:x_mas] = star

	data = fitsio.read(image)
	copy_data = copy.deepcopy(data)
	copy_data_star = copy_data + zeros
	hdu = fits.PrimaryHDU(copy_data_star)
	hdulist = fits.HDUList([hdu])
	hdulist.writeto(out_image_star)

	image_ToO = out_image_star

	os.system("hotpants "\
	 			"-inim {} "\
 				"-tmplim {} "\
 				"-outim {} "\
 				"-c t "\
 				"-tu 50000 "\
 				"-iu 50000 "\
 				"-tl -500 "\
 				"-il -500 "\
 				"-tr 9.5 "\
 				"-ir 9.5".format(image_ToO, image_template, image_diff))

	mag_tup = get_mag_rec(image_diff, CC, circle_radio, x, y) 
	return mag_tup
