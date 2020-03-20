import time
start_time = time.time()

import pandas as pd
import lightkurve as lk
import glob
from astroquery.mast import Catalogs
from scipy import signal
import numpy as np
from astroquery.vizier import Vizier
from astropy import units as u 
import astropy.coordinates as coord
####

def wrms(a, e):
	""" 
	Weighted root mean square of array `a`, with uncertanty given by `e` 
	Parameters
	----------
	a : array
	Array containing data
	e : array
	Uncertainties on `a`.
	The weighted rms is calculated using the weighted mean, where the 
	weights are equal to 1/e**2
	"""
	a = a[np.logical_not(np.isnan(a))]   ## to remove nan in array if any  
	e = e[np.logical_not(np.isnan(e))]
	w = 1/e**2
	return np.sqrt(np.sum(w*(a - np.average(a, weights=w))**2) / sum(w))


def bcflow(teff):
    """Table 1 from Torres 2010: Bolometric Corrections by Flower (1996) as a Function of Temperature:
    BC_V = a + b(log T_eff) + c(log T_eff)^2 + sdot sdot sdot
    Coefficient	log T_eff < 3.70	3.70 < log T_eff < 3.90	   log T_eff>3.90
    a	       -0.190537291496456E+05	-0.370510203809015E+05	-0.118115450538963E+06
    b	        0.155144866764412E+05	 0.385672629965804E+05	 0.137145973583929E+06
    c	       -0.421278819301717E+04	-0.150651486316025E+05	-0.636233812100225E+05
    d	        0.381476328422343E+03	 0.261724637119416E+04	 0.147412923562646E+05
    e	                    ... 	    -0.170623810323864E+03	-0.170587278406872E+04
    f	                    ... 	            ... 	         0.788731721804990E+02
    """

    a = [-0.190537291496456e+05, -0.370510203809015e+05, -0.118115450538963e+06]
    b = [0.155144866764412e+05,   0.385672629965804e+05,  0.137145973583929e+06]
    c = [-0.421278819301717e+04, -0.150651486316025e+05, -0.636233812100225e+05]
    d = [0.381476328422343e+03,   0.261724637119416e+04,  0.147412923562646e+05]
    e = [-0.170623810323864e+03, -0.170587278406872e+04]
    f = [0.788731721804990e+02]

    lteff= np.log10(teff)
    if lteff < 3.7:
        bc = a[0] + (b[0]*lteff) + (c[0]*(lteff**2)) + (d[0]*(lteff**3))
    elif (lteff >= 3.7) and (lteff < 3.9):
        bc = a[1] + (b[1]*lteff) + (c[1]*(lteff**2)) + (d[1]*(lteff**3)) + (e[0]*(lteff**4))
    elif lteff >= 3.9:
        bc = a[2] + (b[2]*lteff) + (c[2]*(lteff**2)) + (d[2]*(lteff**3)) + (e[1]*(lteff**4)) + (f[0]*(lteff)**5)
    return bc


def logg_trigomonetric(teff, mass, v, bc, par, dpar, dteff, dmass):
    """Calculate the trigonometric logg and error"""
    #np.geterr()
    if mass == 'nan':
        logg, dlogg = 'nan', 'nan'
    else:
        e = 2.718281828
        logg  = 4.44 + np.log10(mass) + (4.0*np.log10(teff/5777.)) + (0.4*(v + bc)) + (2.0*np.log10(par/1000.0)) + 0.108
        logg  = np.round(logg, 2)
        dlogg = np.sqrt(((dmass*np.log10(e))/mass)**2 + ((4.*dteff*np.log10(e))/teff)**2 + ((2.*0.05*np.log10(e))/par)**2)
        dlogg = np.round(dlogg, 2)
    return logg, dlogg

# J/AJ/156/102/table9 TESS input catalog and candidate target list
####

# df.query('B > 50 and C != 900')
# data : https://tess.mit.edu/observations/sector-1/



sec1 = pd.read_csv(f'all_targets_S001_v1.txt', sep="\t")
sec2 = pd.read_csv('all_targets_S002_v1.txt', sep="\t")
sec3 = pd.read_csv('all_targets_S003_v1.txt', sep="\t")
sec4 = pd.read_csv(f'all_targets_S004_v1.txt', sep="\t")
sec5 = pd.read_csv('all_targets_S005_v1.txt', sep="\t")
sec6 = pd.read_csv('all_targets_S006_v1.txt', sep="\t")
sec7 = pd.read_csv(f'all_targets_S007_v1.txt', sep="\t")
sec8 = pd.read_csv('all_targets_S008_v1.txt', sep="\t")
sec9 = pd.read_csv('all_targets_S009_v1.txt', sep="\t")
sec10 = pd.read_csv(f'all_targets_S010_v1.txt', sep="\t")
sec11 = pd.read_csv('all_targets_S011_v1.txt', sep="\t")
sec12 = pd.read_csv('all_targets_S012_v1.txt', sep="\t")
sec13 = pd.read_csv('all_targets_S013_v1.txt', sep="\t")

main_cat = pd.concat([sec1,sec2,sec3,sec4,sec5,sec6,sec7,sec8,sec9,sec10,sec11,sec11,sec13]) 

main_cat = main_cat.drop_duplicates(subset='TICID', keep='first') 

main_cat = main_cat.query('Tmag < 9.0')





for i , line in main_cat.iterrows():
	tic_name = 'TIC' + str(int(line['TICID']))
	# catalog_gaia2 = Catalogs.query_object(tic_name, radius=0.001666666666666666, catalog="GaiaDR2")
	# if len(catalog_gaia2) > 0:
	# 	catalog_gaia2 = catalog_gaia2[0]
	# try:
	# 	input_tess_vizier = Vizier.query_region(coord.SkyCoord(ra=catalog_gaia2['ra'], dec= catalog_gaia2['dec'],unit=(u.deg, u.deg), frame='icrs'), width="0.3333336m", catalog=['J/AJ/156/102/table9'])
	# 	mass_tess = input_tess_vizier[0]['M_'][0]
	# 	teff_gaia, teff_err_gaia = catalog_gaia2['teff_val'][0] , abs(catalog_gaia2['teff_percentile_lower'][0] - catalog_gaia2['teff_percentile_upper'][0])
	# 	Vmag_tess = input_tess_vizier[0]['Vmag'][0]
	# 	bc = bcflow(teff_gaia)
	# 	par, dpar = catalog_gaia2['parallax'][0] , catalog_gaia2['parallax_error'][0] 
	# 	log_g = logg_trigomonetric(teff_gaia, mass_tess, Vmag_tess, bc, par, dpar, teff_err_gaia, 0.0)
	# 	print('log g = ', log_g)
	# except :
	# 	log_g = [-1,-1]
	if 'saeed' == 'saeed':
		star_name = line['TICID']
		try:
			t = lk.search_lightcurvefile(star_name).table

			obs_id = t['obs_id']
			
			g = glob.glob(f"/Users/saeedh/.lightkurve-cache/*/*/*/*{obs_id[0][24:40]}*_lc.fits")

			if len(g) > 0:
				print('using cache ... ')
			if len(g) == 0:
				#101010111111
				g = glob.glob(f"/Users/saeedh/.lightkurve-cache/*/*/*/*{obs_id[0][24:40]}*_lc.fits")
				print('downloading data ... ')
			pixels = lk.search_lightcurvefile(star_name).download_all() 
		except:
			print(f'TIC = {star_name} has a problem ... ')

		try: ## TESS
			data = pixels[0]
			lc_pdc = data.PDCSAP_FLUX.remove_nans()
			lc_sap = data.SAP_FLUX.remove_nans().normalize()
			if len(g) == 2:
				for i in np.arange(1,2):
					# data = lk.open(g[i])
					data = pixels[i]
					lc_pdc = lc_pdc.append(data.PDCSAP_FLUX.remove_nans().normalize())
					lc_sap = lc_sap.append(data.SAP_FLUX.remove_nans().normalize())
					print(f'Merging {star} LCs ... {i+1} out of {len(g)}, sector')
			elif len(g) > 2:
				for i in np.arange(1,len(g)):
					# data = lk.open(g[i])
					data = pixels[i]
					lc_pdc = lc_pdc.append(data.PDCSAP_FLUX.remove_nans().normalize())
					lc_sap = lc_sap.append(data.SAP_FLUX.remove_nans().normalize())
					print(f'Merging {star} LCs ... {i+1} out of {len(g)}, sector')
			else: pass

			lc_pdc , lc_sap = lc_pdc.remove_outliers(), lc_sap.remove_outliers(sigma=2) #.flatten(window_length=14177)
			flux_pdc, flux_err_pdc = lc_pdc.flux, lc_pdc.flux_err
			flux_sap_raw, flux_err_sap = lc_sap.flux, lc_sap.flux_err

			flux_sap = signal.savgol_filter(flux_sap_raw,15,1)  # smothing
			print('lc analysing ... ')

			try:
				rms = wrms(flux_sap,flux_err_sap)
			except:
				rms = -1

			print(rms)
			file = open(f'Result_out.rdb', 'a+')
			file.writelines((str(tic_name),',', str(rms*1000), str('\n')))
			file.close()
		except: pass






print(f'FINISH, Time: {time.time() - start_time} sec')