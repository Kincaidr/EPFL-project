
from optparse import OptionParser
import sys
from zwcl2341_plots import zwcl_galaxy_distribution,velocity_dispersion,redshift_plots,radio_SFR_plots
from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import csv

def main( argv ):
    
    SDSS_catalogue = pd.read_csv(r'zwcl2341DR17UncleanViral_kincaid_1.csv')
   
    radio_SDSS_catalogue = pd.read_csv(r'SDSS_radio_combined.csv')
    
    galaxy_ra=SDSS_catalogue['ra']
    galaxy_dec=SDSS_catalogue['dec']
    zsp=SDSS_catalogue['zsp']
    zph=SDSS_catalogue['zph']
    
    radio_SDSS_SFR_catalogue = pd.read_csv(r'solar_mass.csv')

    radio_flux=radio_SDSS_SFR_catalogue['Total_flux']
    solar_mass=radio_SDSS_SFR_catalogue['logMass']
    zsp_1=radio_SDSS_SFR_catalogue['zsp']
    
    cluster_centre= SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))
    
    R_200= 0.103 #5 Mpc

    # zwcl_galaxy_distribution(galaxy_ra,galaxy_dec,zsp,cluster_centre,R_200)
   
    sigma_cluster_z,new_cluster_z=velocity_dispersion(galaxy_ra,galaxy_dec,zsp,cluster_centre,R_200)
    
    redshift_plots(zsp,zph,sigma_cluster_z,new_cluster_z)
    
    #radio_SFR_plots(radio_flux,solar_mass,zsp_1)
    
    
 
     
if __name__ == "__main__":
    main(sys.argv)