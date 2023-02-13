from astroquery.sdss import SDSS
from astropy import coordinates as coords
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.cosmology as cp
import sys
from astropy.io import fits
import pandas as pd
import csv
from astropy.table import Table   
from sklearn import preprocessing
import scipy
from scipy import stats


SDSS = pd.read_csv(r'SDSS_radio_combined.csv')

zph=SDSS['zph']
zsp=SDSS['zsp']

mask_array=~np.logical_or(np.isnan(zsp), np.isnan(zph))


# index_array=np.logical_or(~np.isnan(zph),~np.isnan(zsp))

zph_final=zph[mask_array]
zsp_final=zsp[mask_array]


plt.scatter(zph_final,zsp_final)
plt.xlabel('photometric redshift')
plt.ylabel('spectrscopic redshift')


SDSS_galaxy_ra=np.array(SDSS['ra_1'])
SDSS_galaxy_dec=np.array(SDSS['dec_1'])

radio_galaxy_ra=np.array(SDSS['RA_2'])
radio_galaxy_dec=np.array(SDSS['DEC_2'])

DEC_sep=(SDSS_galaxy_dec-radio_galaxy_dec)*3600

RA_sep=(SDSS_galaxy_ra-radio_galaxy_ra)*3600


#RA and DEC seperation


n, bins = np.histogram(RA_sep, bins=50)

bin_RA_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2


RA_sep_norm = preprocessing.normalize([RA_sep])

n, bins_RA,_ = plt.hist(RA_sep_norm,bins=50,density=1)

mu, sigma = scipy.stats.norm.fit(RA_sep_norm)
best_fit_line_RA = scipy.stats.norm.pdf(bins_RA, mu, sigma)

plt.plot(bins_RA, best_fit_line_RA)

plt.show()

fig1 = plt.figure()
ay1 = fig1.add_subplot(1,1,1, aspect='equal')
ay1.set_xlabel("RA seperation (arcsec)")
ay1.set_title("RA seperation plot for radio and optical sources")
plt.plot(bins_RA, best_fit_line_RA)

#n, bins, patches = plt.hist(DEC_sep, bins=50, color='black', alpha=0.1, rwidth=0.85)

n, bins = np.histogram(DEC_sep, bins=50,density=1)

bin_DEC_centre=(bins[np.argmax(n)]+bins[np.argmax(n)+1])/2

n, bins_DEC = np.histogram(DEC_sep-bin_DEC_centre,bins=50,density=1)

mu, sigma = scipy.stats.norm.fit(DEC_sep-bin_DEC_centre)
best_fit_line_DEC = scipy.stats.norm.pdf(bins_DEC, mu, sigma)

fig2 = plt.figure()
ay2 = fig2.add_subplot(1,1,1, aspect='equal')
ay2.set_xlabel("DEC seperation (arcsec)")
ay2.set_title("DEC seperation plot for radio and optical sources")
plt.plot(bins_DEC,best_fit_line_DEC)
plt.hist(bins)
plt.show()


        




# input_table = Table.read('cross_match.csv')

# table = XMatch.query(cat1=input_table,cat2='SIMBAD',max_distance=5 * u.arcsec, colRA1='ra',colDec1='dec')


# table = XMatch.query(cat1='zwcl_radio.csv',cat2='vizier:II/246/out',
# max_distance=5 * u.arcsec, colRA1='RA',colDec1='DEC',colRA2='RA_ICRS',colDec2='DE_ICRS')


# print(table['main_type'])

# print('length of cross match:',len(cross_match_ra))
# header=['ra','dec']
# data=[cross_match_ra,cross_match_dec]

# with open('cross_match.csv', 'w', encoding='UTF8', newline='') as f:
#     writer = csv.writer(f)
#     writer.writerow(header)
#     for i in range(len(cross_match_ra)):
#         data=[cross_match_ra[i],cross_match_dec[i]]
#         writer.writerow(data)
        
# print('arcmin circle',len(arcmin_circle))
#  if np.sqrt((galaxy_ra[i]-cluster_centre[0])**2+(galaxy_dec[i]-cluster_centre[1])**2)<= (Angle('5arcmin',u.degree).value):
# cross_match=[]
# for i in range(len(galaxy_ra_radio)):
#     for j in range(len(galaxy_ra_SDSS)):
#         if np.sqrt((galaxy_ra_radio[i]-galaxy_ra_SDSS[j])**2+(galaxy_dec_radio[i]-galaxy_dec_SDSS[j])**2)<= (Angle('1arcmin',u.degree).value):
#             cross_match.append(galaxy_ra_radio)

# galaxy_ra_SDSS=np.array(SDSS['ra'])
# galaxy_dec_SDSS=np.array(SDSS['dec'])

# galaxy_ra_radio=np.array(radio['RA'])
# galaxy_dec_radio=np.array(radio['DEC'])

# # cross match

# SDSS_catalogue=SkyCoord(Angle(galaxy_ra_SDSS,u.degree),Angle(galaxy_dec_SDSS,u.degree))
# radio_catalogue=SkyCoord(Angle(galaxy_ra_radio,u.degree),Angle(galaxy_dec_radio,u.degree))


# idx_radio, idx_SDSS, d2d, d3d = SDSS_catalogue.search_around_sky(radio_catalogue, 4*u.arcsec)

# cross_match_ra=SDSS_catalogue[idx_SDSS].ra.value

# cross_match_dec=SDSS_catalogue[idx_SDSS].dec.value



# DEC_sep=(SDSS_catalogue[idx_SDSS].dec.value-radio_catalogue[idx_radio].dec.value)*3600

# RA_sep=(SDSS_catalogue[idx_SDSS].ra.value-radio_catalogue[idx_radio].ra.value)*3600