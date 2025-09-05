from rvspecfit import fitter_ccf, vel_fit, spec_fit, utils
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import os

config=utils.read_config('/Users/mncavieres/Documents/2024-1/Delve/xshooter/config.yaml') # optional
# you can create a configuration file with various options


run_path = '/Users/mncavieres/Documents/2024-1/Delve/xshooter/Spec/Star01'
# let's assume we have data stored in 1d arrays in a table
# with wavelength being in Angstrom
# espec being a vector of standard deviations

tab=Table().read('/Users/mncavieres/Documents/2024-1/Delve/xshooter/Spec/Star01/ADP.2025-07-13T14_06_10.099.fits', hdu=1) # read the fits file
wavelength = tab['WAVE']
spec = tab['FLUX_REDUCED'] #+ 2
espec = tab['ERR_REDUCED'] # add a small value to avoid division by zero

# reduce range
lim_mask = (wavelength.value > 1200) & (wavelength.value < 2478)
wavelength = wavelength[lim_mask]
spec = spec[lim_mask]
espec = espec[lim_mask]

espec = np.where(espec <= 0, 1e10, espec)


# This constructs the specData object from wavelength, spectrum and error
# spectrum arrays. The rvspecfit works on arrays of SpecData's which may
# represent multiple exposures or multiple spectral configurations
# Here we just have one spectrum from the spectral configuration "mysetup"

specdata = [spec_fit.SpecData('xshooter',
                               wavelength*10, # convert nm to Angstrom
                               spec,
                               espec#,
                               #badmask=badmask
                               )
                              ]


# this tries to get a sensible guess for the stellar parameters/velocity
res = fitter_ccf.fit(specdata, config)
paramDict0 = res['best_par']
print(paramDict0)

print('Initial guess for CCF: ', res['best_vel'], res['best_vsini'], res['best_par'])


fixParam = [] 
if res['best_vsini'] is not None:
    paramDict0['vsini'] = res['best_vsini']

# start off with a higher feh
paramDict0['feh'] = 0.5


options = {'npoly':15}

# this does the actual fitting performing the maximum likelihood fitting of the data
res1 = vel_fit.process(specdata,
                           paramDict0,
                           fixParam=fixParam,
                           config=config,
                           options=options)
print(res1)

# save the results to a pickle file
import pickle

with open("/Users/mncavieres/Documents/2024-1/Delve/xshooter/Spec/Star01/res_vacuum_default_full.pkl", "wb") as f:
    pickle.dump(res1, f)


# change the badmask on the Halpha region to true
#badmask[(wavelength.value > 6540) & (wavelength.value < 6620)] = False

# # same plot but with matplotlib
plt.figure(figsize=(30/1.5, 10/1.5))
#plt.plot(wavelength[~badmask], spec[~badmask], color='black', label='Spectrum')
plt.plot(wavelength, spec, color='black', label='Spectrum')
plt.plot(wavelength, res1['yfit'][0], color='red', label='Fit')
plt.xlabel('Wavelength [nm]')
plt.ylabel('flux')
plt.title('1D Coadded Spectrum with Fit')
plt.legend()
plt.grid()
plt.savefig(os.path.join(run_path,'rvspecfit_results_NIR.pdf'), bbox_inches='tight')
plt.show()
#Table(res1).write('b576_rvspecfit_results.fits', overwrite=True)
