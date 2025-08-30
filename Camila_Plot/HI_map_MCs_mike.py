import pandas as pd
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
import numpy as np
import gala
from gala.coordinates import MagellanicStreamNidever08
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from matplotlib import patheffects
from astropy.coordinates import SkyCoord
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

df = pd.read_csv('/Users/mncavieres/Documents/2024-1/Delve/Camila_Plot/MagE_targets_x_GaiaBPRP_RVextended.csv')
#df_all = Table.read('/Users/mncavieres/Documents/2024-1/Delve/Proposal/Fase2_Catalog/fase2_catalog.fits').to_pandas()
giants = Table.read('/Users/mncavieres/Documents/2024-1/Delve/Proposal/mike/filtered_stars_with_airmass2.csv')
ra_vc, dec_vc, vgsr_vc = np.genfromtxt('/Users/mncavieres/Documents/2024-1/Delve/Camila_Plot/ChandraMS.txt', usecols=(1,2,4), unpack=True)

# select only stars with min_airmass between 1.0 and 2
giants = giants[(giants['min_airmass'] < 2.0)]

c = coord.FK5(ra=df['RA2000'].values*u.deg, dec=df['DEC2000'].values*u.deg)
ms = c.transform_to(MagellanicStreamNidever08())

#c_all = coord.FK5(ra=df_all['RA2000'].values*u.deg, dec=df_all['DEC2000'].values*u.deg)
#ms_all = c_all.transform_to(MagellanicStreamNidever08())

c_vc = coord.FK5(ra=ra_vc*u.deg, dec=dec_vc*u.deg)
ms_vc = c_vc.transform_to(MagellanicStreamNidever08())

c_lmc = coord.FK5(ra=(15.*(5.+40./60.+5./3600.))*u.deg, dec=(-69.-45./60.-51./3600.)*u.deg)
ms_lmc = c_lmc.transform_to(MagellanicStreamNidever08())

c_smc = coord.FK5(ra=(15.*(0.+52./60.+44.8/3600.))*u.deg, dec=(-72.-49./60.-43./3600.)*u.deg)
ms_smc = c_smc.transform_to(MagellanicStreamNidever08())

c_giants = coord.FK5(ra=giants['RA2000']*u.deg, dec=giants['DEC2000']*u.deg)
ms_giants = c_giants.transform_to(MagellanicStreamNidever08())

fits_file = '/Users/mncavieres/Documents/2024-1/Delve/Camila_Plot/Stream_clean (1).fits'
hdul = fits.open(fits_file)

wcs = WCS(hdul[1].header)

image_data = hdul[1].data
hdul.close()

plt.figure(figsize=(10, 5), frameon=False)
plt.subplot(111, projection="aitoff")

frame = "icrs"

ax = plt.subplot(projection=wcs, frame_class=EllipticalFrame, label='HI map')
im = ax.imshow(image_data, cmap='Greys', origin='lower')

cbr_t= plt.scatter(giants['RA2000'], giants['DEC2000'], transform=ax.get_transform('icrs'),
           s=5, marker='.', c = giants['min_airmass'], label=r'Giants (60 - 100 kpc)')

plt.scatter(ra_vc, dec_vc, transform=ax.get_transform('icrs'),
           s=50, edgecolor='red', facecolor='none', label='Chandra+23')
path_effects=[patheffects.withStroke(linewidth=3, foreground='black')]

#plt.scatter(df_all['RA2000'].values, df_all['DEC2000'].values, transform=ax.get_transform('icrs'),
#           s=50, edgecolor='blue', facecolor='none', label='Targets')
ax.coords['glat'].set_ticklabel(color='w', fontsize=20, zorder=100)
ax.coords.grid(True, color='k', ls='dotted')  # Add a grid for better reference
plt.legend(loc=1, fontsize=14)
ax.coords['glat'].set_ticklabel(color='k', fontsize=20, zorder=100)
plt.colorbar(cbr_t, label='Minimum Airmass')
plt.tight_layout()
plt.savefig('/Users/mncavieres/Documents/2024-1/Delve/Plots/mike/test.pdf', format='pdf')
plt.close()