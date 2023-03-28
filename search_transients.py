import fitsio
import search_transients_utils as stu
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii

frb_name = ["frb180924","frb180924","frb181112","frb181112","frb190102","frb190102","frb190608","frb190608","frb190608",'frb190608',"frb190611","frb190611",'frb190611',"frb190611","frb190711",'frb190711', "frb190711", "frb190711", "frb190714", "frb190714",'frb190714', "frb190714", "frb191001", "frb191001", "frb191001", "frb191001",'frb191001','frb191001']
date_ToO = ['190531','200629','190531','200629','190531','200629','190608_02','190808','190808','200629','190611','190808','190808','200629','190711','190808','190808','200629','190714','190808','190808','200630','191003','191126','191126','200629']
date_template = ['200629','190531','200629','190531','200629','190531','200629','200629','190608_02','190608_02','200629','200629','190611','190611','200629','200629','190711','190711','200630','200630','190714','190714','200629','200629','191003','191003']
FWHM = [5.17, 4.33, 3.92, 9.21, 5.95, 5.28, 5.87, 5.7, 5.7, 8.62, 4.91, 6.17, 6.17, 8.96, 6.06, 6.13, 6.13, 4.87, 4.37, 9.08, 9.08, 4.24, 5.43, 5.9, 5.9, 8.83]

# frb_name = ["frb191001"]
# date_ToO = ['191003']
# date_template = ['191001']
# FWHM = [5.43]


main_tab = Table()
main_tab['frb_name'] = frb_name
main_tab['date_ToO'] = date_ToO
main_tab['date_template'] = date_template
main_tab['FWHM'] = FWHM
x = 2048
y = 2048

upper_lims = []

for frb_name, date_ToO, date_template, FWHM in zip (frb_name, date_ToO, date_template, FWHM):

	image = "/home/consuelo/projects/FRBs/LCOgt/coadds_frb_center/" + frb_name + "/" + frb_name + "_" + date_ToO + "_" + "rp.fits"
	skymapper = '/home/consuelo/projects/FRBs/LCOgt/photometry/mag_galaxies_host/catalogos/' + frb_name + '_skymapper.tsv'
	out_image_star = r'/home/consuelo/projects/FRBs/LCOgt/upper_lims/' + frb_name + "_" + date_ToO + '_rp_upperlimstar_v2.fits'
	filtro = 'r_petro'
	h_size = 10
	#FWHM = 6
	circle_radio = 2 * FWHM
	image_template = "/home/consuelo/projects/FRBs/LCOgt/coadds_frb_center/" + frb_name + "/" + frb_name + "_" + date_template + "_rp.fits"
	image_diff = "/home/consuelo/projects/FRBs/LCOgt/upper_lims/diff_" + frb_name + "_" + date_ToO + "_" + date_template + "_rp_upperlimstar_v2.fits"

	mag_array = np.arange(19,25,0.1)
	mag_diff = []
	err_mag = []
	mag_recov = []
	SNoise = []
	for mag in mag_array:
		mag_rec, err_mag_p, err_mag_m, SN = stu.mag_rec(image, skymapper, out_image_star, filtro, mag, x, y, h_size, FWHM, circle_radio, image_template, image_diff)
		mag_recov += [mag_rec]
		SNoise += [SN]
		err_mag += [np.mean([err_mag_p, err_mag_m])]
		mag_diff += [mag_rec - mag]
		os.remove(out_image_star)
		os.remove(image_diff)
		print('SN = '+str(SN))
		print('mag_add = '+str(mag))
		print('mag_recov = '+ str(mag_rec))

	tab = Table()
	tab['mag'] = mag_array
	tab['mag_rec'] = mag_recov
	tab['mag_diff'] = mag_diff
	tab['SN'] = SNoise

	# SN mayor a 3 => detectada, upper limit maximo de los SN <= 3
	SNoise = np.array(SNoise)
	cond = np.where(SNoise < 3.)
	i = np.min(cond)                   # minimo porque np.where me da los indices
	upper_lim = tab['mag'][i]
	upper_lims += [upper_lim]

main_tab['upper_lim'] = upper_lims


# plt.errorbar(mag_array, mag_diff, yerr = err_mag, fmt = 'ro')
# plt.xlabel('magnitude')
# plt.ylabel('recovered magnitude - magnitude')
# plt.title('Magnitude in position x = {}, y = {}'.format(x, y))
# plt.show()


# # guarda imagen de diferencia con estrella artificial de magnitud limite superior
# stu.mag_rec(image, skymapper, out_image_star, filtro, upper_lim, x, y, h_size, FWHM, circle_radio, image_template, image_diff)


# print(tab)
# ascii.write(tab, '/home/consuelo/Desktop/up/' + frb_name + '_' + date_template + '.dat')
# print('upper limit: {}'.format(upper_lim))