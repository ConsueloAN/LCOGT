import os

frb_name = 'frb191001'
date_ToO = '200629'
date_template = '191003'

image_ToO = "/home/consuelo/projects/FRBs/LCOgt/coadds_frb_center/" + frb_name + "/" + frb_name + "_" + date_ToO + "_rp.fits"
image_template = "/home/consuelo/projects/FRBs/LCOgt/coadds_frb_center/" + frb_name + "/" + frb_name + "_" + date_template + "_rp.fits"
image_diff = "/home/consuelo/projects/FRBs/LCOgt/differences/diff_" + frb_name + "_" + date_ToO + "_" + date_template + "_rp.fits"

os.system("hotpants "\
 			"-inim {} "\
			"-tmplim {} "\
			"-outim {} "\
			"-tu 50000 "\
			"-iu 50000 "\
			"-tl -500 "\
			"-il -500 "\
			"-tr 9.5 "\
			"-ir 9.5".format(image_ToO, image_template, image_diff))
