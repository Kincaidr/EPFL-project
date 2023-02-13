from astroquery.sdss import SDSS
from astropy import coordinates as coords
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
import numpy as np
from astropy.stats import biweight_scale
from astropy.io import ascii
from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = -1
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.pyplot import show, plot, draw
from multiprocessing import Process
from scipy.constants import pi
import scipy.constants as constants
import statistics as stats
import astropy.cosmology as cp
import sys
from astropy.io import fits
from matplotlib.patches import Circle
import pandas as pd
import csv

def histogram(sample,bins,color):
  n, bins, patches = plt.hist(sample, bins=bins, color=color, alpha=0.1, rwidth=0.85)
  cluster_z=bins[np.argmax(n)]
  return cluster_z

pos = coords.SkyCoord('23h43m39.700s +0d19m51.000s', frame='icrs')
xid = SDSS.query_region(pos, radius='1deg',spectro=True,data_release=17,photoobj_fields=['objid','ra','dec','u','g','r','i'],specobj_fields=['z','class'])


query='SELECT TOP 10 * FROM PhotoObj AS p, stellarMassFSPSGranEarlyDust AS sf, SpecObj AS s JOIN dbo.fGetNearbyObjEq(355.91541666667,0.33083333333333,60) AS A ON (s.bestobjid = A.objid)  WHERE (s.z > 0 AND s.z < 1) AND (s.class= "GALAXY" OR s.class="QSO") AND (s.instrument = "SDSS" OR s.instrument = "BOSS")'

xid = SDSS.query_sql(query)


df=xid.to_pandas()
df.to_csv('./zwcl_SDSS.csv')

#xid = pd.read_csv(r'zwcl2341_90_DR17.csv')

#xid.colnames

redshift=[]
galaxy_ra=[]
galaxy_dec=[]


for i in range(0,len(xid['z'])):
    
    if ((xid['z'][i]<=0.29) and (xid['z'][i]>=0.26)):
        
        redshift.append(xid['z'][i])
        galaxy_ra.append(xid['ra'][i])
        galaxy_dec.append(xid['dec'][i])
        
        
header=['ra','dec','z']
data=[galaxy_ra,galaxy_dec,redshift]

with open('zwcl_SDSS.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    for i in range(len(galaxy_ra)):
        data=[galaxy_ra[i],galaxy_dec[i],redshift[i]]
        writer.writerow(data)
        
redshift_filter=np.array(redshift)

SDSS = pd.read_csv(r'zwcl2341_45_DR17.csv')

csv_reader = csv.reader(SDSS, delimiter = ',')

list_of_column_names=[]

for row in csv_reader:
  list_of_column_names.append(row)


cluster_z=histogram(redshift,500,'#0504aa')
cluster_z_n=histogram(redshift,100,'#0504aa')
cluster_z_median= np.median(redshift)

print('cluster z 500 bin size:',cluster_z)
print('cluster z 100 bin size:',cluster_z_n)
print('cluster z median:',cluster_z_median)

print('total number of galaxies in SDSS is',len(redshift))
print('min redshift in SDSS is',min(redshift))
print('max redshift in SDSS is',max(redshift))            
print('RA and DEC of galaxies in redshift range:',len(galaxy_ra))



#n, bins, patches = plt.hist(redshift, bins=30, color='#0504ab', alpha=0.1, rwidth=0.85)

# cluster_z=bins[np.argmax(n)]
# print('cluster z large sample:',cluster_z)

 
ra_radian = np.array(galaxy_ra)*np.pi/180

dec_radian = np.array(galaxy_dec)*np.pi/180


h=0.7
H0=h*100
cosmo = cp.FlatLambdaCDM(H0=h*100, Om0=0.30)


sc = SkyCoord(ra_radian, dec_radian, unit='rad', representation_type='unitspherical')
cartesian=sc.cartesian


x_coord=(np.cos(dec_radian) * np.cos(ra_radian))
y_coord=(np.cos(dec_radian) * np.sin(ra_radian))
z_coord=(np.sin(dec_radian))

cluster_centre = [355.91541666667,0.33083333333333]

# x_coord=np.array(cartesian.x)
# y_coord=np.array(cartesian.y)
# z_coord=np.array(cartesian.z)

z_filter=np.array(redshift_filter)


cluster_centre= SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))

cluster_centre_ra_radian=cluster_centre.ra.radian
cluster_centre_dec_radian=cluster_centre.dec.radian

cluster_centre=[cluster_centre.ra.value,cluster_centre.dec.value]


sc_cluster_centre=  SkyCoord(cluster_centre_ra_radian, cluster_centre_dec_radian, unit='rad', representation_type='unitspherical')
cartesian_cluster_centre=sc_cluster_centre.cartesian

cluster_centre_x=cartesian_cluster_centre.x.value
cluster_centre_y=cartesian_cluster_centre.y.value
cluster_centre_z=cartesian_cluster_centre.z.value

comoving_centre_x=cosmo.comoving_distance((cluster_z_median)*cluster_centre_x).value
comoving_centre_y=cosmo.comoving_distance((cluster_z_median)*cluster_centre_y).value
comoving_centre_z=cosmo.comoving_distance((cluster_z_median)*cluster_centre_z).value


# reference to centre

comoving_x=((cosmo.comoving_distance(z_filter)*x_coord).value)
comoving_y=((cosmo.comoving_distance(z_filter)*y_coord).value)
comoving_z =((cosmo.comoving_distance(z_filter)*z_coord).value)



R_200= 5


# R_200= 0.1388
      
fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
ay1 = fig1.add_subplot(1,1,1, aspect='auto',projection='3d')
ay2 = fig2.add_subplot(1,1,1, aspect='equal')
ay3 = fig3.add_subplot(1,1,1, aspect='equal')

# 2D physical coordinates

ay1.scatter(comoving_x, comoving_y, comoving_z, color='black', alpha=0.5)
ay1.set_xlabel("X in Mpc")
ay1.set_ylabel("Y in Mpc")
ay1.set_zlabel("Z in Mpc")
ay1.set_title("3-D distribution of galaxies in a 1.5 deg radius centered on Zwcl 2341")

ay2.scatter(comoving_y, comoving_z,color='black', alpha=0.5)
ay2.set_xlabel("Y in Mpc")
ay2.set_ylabel("Z in Mpc")


circle=plt.Circle((comoving_centre_y, comoving_centre_z), R_200, edgecolor= 'blue',
facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed') #, )
ay2.add_patch(circle)
ay2.text(0,0, "R_200", ha="left", va="top",fontsize=10)
fig2.savefig("Y_Z_Mpc")
print(comoving_centre_y,comoving_centre_z)
R_200= 0.103

# 2D angular coordinates

ay3.scatter(galaxy_ra, galaxy_dec,color='black', alpha=0.5)

ay3.set_xlim(min(galaxy_ra),max(galaxy_ra))
ay3.set_ylim(min(galaxy_dec),max(galaxy_dec))
ay3.set_xlabel("RA in deg")
ay3.set_ylabel("DEC in deg")
ay3.set_title("1.5 deg SDSS DR17 galaxies of Zwcl 2341 (Shishir catalogue)")

circle=plt.Circle((cluster_centre[0], cluster_centre[1]), R_200, edgecolor= 'blue',
facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed') #, )
ay3.add_patch(circle)
ay3.text(0,0, "R_200", ha="left", va="top",fontsize=10)
fig3.savefig("Y_Z_Mpc")


ra_5Mpc_circle=[]
dec_5Mpc_circle=[]
z_5Mpc_circle=[]

for i in range(len(galaxy_ra)):
    if np.sqrt((galaxy_ra[i]-cluster_centre[0])**2+(galaxy_dec[i]-cluster_centre[1])**2)<= R_200:
        ra_5Mpc_circle.append(galaxy_ra[i])
        dec_5Mpc_circle.append(galaxy_dec[i])
        z_5Mpc_circle.append(redshift[i])

print('arcmin circle',len(z_5Mpc_circle))

  
new_cluster_z=histogram(z_5Mpc_circle,50,'#0504ab')


print('new cluster z:' ,new_cluster_z)  



cluster_z_filter=[]

for i in range(len(z_5Mpc_circle)):
    if (z_5Mpc_circle[i] >= new_cluster_z - 0.015) and (z_5Mpc_circle[i] <= new_cluster_z + 0.015):
        cluster_z_filter.append(z_5Mpc_circle[i]) 
print(max(cluster_z_filter),min(cluster_z_filter),len(cluster_z_filter))


  
  #Cluster members
  
 
v_rec_vel=[]         
for i in range(len(cluster_z_filter)):
   v_rec_vel.append(cluster_z_filter[i]*constants.c)
    
median_rec_vel=stats.median(v_rec_vel)
        
v_rest_vel=[]
for i in range(len(cluster_z_filter)):
    v_rest_vel.append(((v_rec_vel[i] - median_rec_vel) /(1+median_rec_vel/constants.c)) )

         
biscl_2 = biweight_scale(v_rest_vel)

rest_vel_sigma_clip=[]

for i in range(len(v_rest_vel)):
  if ((-3*biscl_2 ) <= v_rest_vel[i]) and (3*biscl_2  >= v_rest_vel[i]):
    rest_vel_sigma_clip.append(v_rec_vel[i])

print(' rest velocity calculated from the biweight scale',(biscl_2/10**3))

biscl_2_sigm_clip = biweight_scale(rest_vel_sigma_clip)

print(' rest velocity calculated from the biweight scale sigma clip',(biscl_2_sigm_clip/10**3))

plt.show() 

