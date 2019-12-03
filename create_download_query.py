
# Create a download .sh file of TESS LC for a list of TESS target based on TIC id
#  
#


# https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html

from astroquery.mast import Observations
import pandas as pd
import numpy as np
import glob

#####

list_cat = glob.glob("./store/all_target*.txt") # read the list of all TIC id within sector 1:13

for i in np.arange(len(list_cat)):
	if i < 1:
		main_cat = pd.read_csv(list_cat[0], sep="\t")
	else:
		sec = pd.read_csv(list_cat[i], sep="\t")
		main_cat = main_cat.append(sec)

main_cat = main_cat.drop_duplicates(subset='TICID', keep='first') # remove repeted target in the different sectors

main_cat = main_cat.query('Tmag < 9.0') # this is a list of all TIC id in sectors 1:13, filtered by Tmag < 9.0

for i in np.arange(1,14):
	down_cat = pd.read_csv(f'./store/tesscurl_sector_{i}_lc.sh', delim_whitespace=False, header=None, skiprows=[0])
	for i , tar in main_cat.iterrows():
		down_cat['index'] = down_cat[0].str.find(str(int(tar['TICID'])))
		s = down_cat[down_cat['index'] != -1][0]
		if len(s) > 0:
			s.to_csv('./Data_storage/final_download.sh', mode='a', index=None)

print('Done !')




