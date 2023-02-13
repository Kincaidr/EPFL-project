
import astropy.io.fits as fits
import requests 


url='https://www.legacysurvey.org/viewer/cutout.fits?ra=355.9186&dec=0.2747&layer=hsc-dr2&pixscale=5.00&size=1100' 
r=requests.get(url) 


fits_outname='test.fits'

open(fits_outname,"wb").write(r.content)

data = fits.open(fits_outname)[0].data




