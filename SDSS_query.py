from astroquery.sdss import SDSS
import pandas as pd
from astropy import coordinates as coords
import csv

def histogram(sample,bins,color):
  n, bins, patches = plt.hist(sample, bins=bins, color=color, alpha=0.1, rwidth=0.85)
  cluster_z=bins[np.argmax(n)]
  return cluster_z

pos = coords.SkyCoord('23h43m39.700s +0d19m51.000s', frame='icrs')
xid = SDSS.query_region(pos, radius='1deg',spectro=True,data_release=17,photoobj_fields=['objid','ra','dec','u','g','r','i'],specobj_fields=['z','class'])


query_galaxy_SDSS='SELECT TOP 10000 * FROM PhotoObj AS p, stellarMassFSPSGranEarlyDust AS sf, SpecObj AS s JOIN dbo.fGetNearbyObjEq(355.91541666667,0.33083333333333,60) AS A ON (s.bestobjid = A.objid)  WHERE (s.z > 0 AND s.z < 1) AND s.class= "GALAXY"  AND s.instrument = "SDSS" '
query_galaxy_BOSS='SELECT TOP 10000 * FROM PhotoObj AS p, stellarMassFSPSGranEarlyDust AS sf, SpecObj AS s JOIN dbo.fGetNearbyObjEq(355.91541666667,0.33083333333333,60) AS A ON (s.bestobjid = A.objid)  WHERE (s.z > 0 AND s.z < 1) AND s.class= "GALAXY" AND  s.instrument = "BOSS"'
query_QSO_SDSS='SELECT TOP 10000 * FROM PhotoObj AS p, stellarMassFSPSGranEarlyDust AS sf, SpecObj AS s JOIN dbo.fGetNearbyObjEq(355.91541666667,0.33083333333333,60) AS A ON (s.bestobjid = A.objid)  WHERE (s.z > 0 AND s.z < 1) AND s.class="GALAXY" AND s.instrument = "SDSS"'
query_QSO_BOSS='SELECT TOP 10000 * FROM PhotoObj AS p, stellarMassFSPSGranEarlyDust AS sf, SpecObj AS s JOIN dbo.fGetNearbyObjEq(355.91541666667,0.33083333333333,60) AS A ON (s.bestobjid = A.objid)  WHERE (s.z > 0 AND s.z < 1) AND  s.class="QSO" AND s.instrument = "BOSS"'
xid_galaxy_SDSS = SDSS.query_sql(query_galaxy_SDSS,timeout=1200)
xid_galaxy_BOSS = SDSS.query_sql(query_galaxy_BOSS,timeout=1200)
xid_QSO_SDSS = SDSS.query_sql(query_QSO_SDSS,timeout=1200)
xid_QSO_BOSS = SDSS.query_sql(query_QSO_BOSS,timeout=1200)

df=xid_galaxy_SDSS.to_pandas()
df.to_csv('./zwcl_galaxy_SDSS.csv')

df=xid_galaxy_SDSS.to_pandas()
df.to_csv('./zwcl_galaxy_BOSS.csv')

df=xid_QSO_SDSS.to_pandas()
df.to_csv('./zwcl_QSO_SDSS.csv')

df=xid_QSO_BOSS.to_pandas()
df.to_csv('./zwcl_QSO_BOSS.csv')
