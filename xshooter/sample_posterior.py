from rvspecfit import fitter_ccf, vel_fit, spec_fit, utils
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
import os
from multiprocessing import Pool


plot_dir = '/Users/mncavieres/Documents/2024-1/Delve/xshooter/Spec/Star01/VIS_FIT'

#config=utils.read_config('/Users/mncavieres/Documents/2025-1/B576/rvspecfit_local/config.yaml') # optional
# you can create a configuration file with various options

# helper functions
def save_traces(sampler, name):
    chain = sampler.get_chain()
    fig, axes = plt.subplots(ndim,1,figsize=(10,2.5*ndim), sharex=True)
    for i in range(ndim):
        for w in range(nwalkers):
            axes[i].plot(chain[:,w,i], alpha=0.4, lw=0.7)
        axes[i].set_ylabel(labels[i])
        axes[i].grid(True)
    axes[-1].set_xlabel('Step')
    plt.tight_layout()
    fig.savefig(os.path.join(plot_dir, f'trace_{name}.pdf'))
    plt.close(fig)


# let's assume we have data stored in 1d arrays in a table
# with wavelength being in Angstrom
# espec being a vector of standard deviations
config=utils.read_config('/Users/mncavieres/Documents/2024-1/Delve/xshooter/config.yaml') # optional
# you can create a configuration file with various options


run_path = '/Users/mncavieres/Documents/2024-1/Delve/xshooter/Spec/Star01'
# let's assume we have data stored in 1d arrays in a table
# with wavelength being in Angstrom
# espec being a vector of standard deviations

tab=Table().read('/Users/mncavieres/Documents/2024-1/Delve/xshooter/Spec/Star01/ADP.2025-07-13T14_06_10.094.fits', hdu=1) # read the fits file
wavelength = tab['WAVE']
spec = tab['FLUX_REDUCED'] #+ 2
espec = tab['ERR_REDUCED'] # add a small value to avoid division by zero

# reduce range
lim_mask = (wavelength.value > 552) & (wavelength.value < 1019)
wavelength = wavelength[lim_mask]
spec = spec[lim_mask]
espec = espec[lim_mask]

#mask anything with espec <=0, by settin espec to a huge value

espec = np.where(espec <= 0, 1e10, espec)

# also mask between 757.5 and 767.5 nm (O2 telluric)
espec = np.where((wavelength.value > 757.5) & (wavelength.value < 767.5), 1e10, espec)

# This constructs the specData object from wavelength, spectrum and error
# spectrum arrays. The rvspecfit works on arrays of SpecData's which may
# represent multiple exposures or multiple spectral configurations
# Here we just have one spectrum from the spectral configuration "mysetup"

specdata = [spec_fit.SpecData('vis_xshooter',
                               wavelength*10, # convert nm to Angstrom
                               spec,
                               espec#,
                               #badmask=badmask
                               )
                              ]



# This constructs the specData object from wavelength, spectrum and error
# spectrum arrays. The rvspecfit works on arrays of SpecData's which may
# represent multiple exposures or multiple spectral configurations
# Here we just have one spectrum from the spectral configuration "mysetup"

# specdata = [spec_fit.SpecData('vus_xshooter',
#                             wavelength*10, # nm to Angstrom
#                             spec,
#                             espec#,
#                             #badmask=badmask
#                             )
#                             ]


# this tries to get a sensible guess for the stellar parameters/velocity
res = fitter_ccf.fit(specdata, config)
paramDict0 = res['best_par']
print(paramDict0)
# set the initial guess for the parameters using Kreuzer+2020
#paramDict0['teff'] = 11400 #res['best_teff']
#paramDict0['logg'] = 3.79
#paramDict0['feh'] = 0 #res['best_feh']


fixParam = [] 
if res['best_vsini'] is not None:
    paramDict0['vsini'] = res['best_vsini'] # 47 as per Kreuzer+2020

options = {'npoly':15}

# this does the actual fitting performing the maximum likelihood fitting of the data
res1 = vel_fit.process(specdata,
                        paramDict0,
                        #fixParam=fixParam,
                        config=config,
                        options=options)

# with this inital fit we get the radial velocity now lets use the likelihood to sample the posterior distribution
# for the Teff, logg, feh and alpha parameters

#loglike = 0
vel = res1['vel']
print('The best fit radial velocity is: ', vel)


ndim, nwalkers, nsteps = 4, 400, 3000#4, 400, 3000

# initial guess for the parameters, flat prior
teff_init = np.random.uniform(2000, 15000, nwalkers)
logg_init = np.random.uniform(0.1, 6.2, nwalkers)
feh_init = np.random.uniform(-4,1, nwalkers)
alpha_init = np.random.uniform(-0.19, 1.1, nwalkers)

p0 = np.array([teff_init, logg_init, feh_init, alpha_init]).T



# define the log likelihood function
def log_likelihood(theta, specdata, vel):
    teff, logg, feh, alpha = theta

    #vel=300
    atm_params = [teff,logg, feh, alpha] # 'teff', 'logg', 'feh', 'alpha'
    spec_chisq  = spec_fit.get_chisq(specdata,
                                vel,
                                atm_params,
                                None,
                                config=config, options = dict(npoly=15))

    return -0.5 * spec_chisq # return the negative log likelihood

if __name__ == "__main__":


    # set up pool for parallel processing
    with Pool(8) as pool:
        sam_spec = emcee.EnsembleSampler(
            nwalkers, ndim, log_likelihood,
            args=(specdata, vel), vectorize=False, pool=pool
            )

        sam_spec.run_mcmc(p0, nsteps, progress=True)


    flat_spec = sam_spec.get_chain(discard=1000, flat=True)
    labels = ['Teff','logg','[Fe/H]','alpha']

    np.save(os.path.join(plot_dir, 'samples_spec.npy'), flat_spec)

    # corner plots courtesy of chatgpt
    fig_c = corner.corner(
        flat_spec, labels=labels, color='k', fill_contours=False,
        levels=[0.68,0.95,0.997], contour_kwargs={'linewidths':1.5}, contourf_kwargs={'alpha':0.3}, smooth=1, bins=30, 
        show_titles=True, title_kwargs={"fontsize": 12},
    )

    # fig_c.text(0.65, 0.92, 'Green: Photometry-only', color='green', fontsize=12)
    # fig_c.text(0.65, 0.88, 'Blue: Photometry+Spectroscopy', color='blue', fontsize=12)
    plt.tight_layout()
    fig_c.savefig(os.path.join(plot_dir,'corner_spec.pdf'))
    plt.close(fig_c)

    # plot the traces

    # trace plots
    save_traces(sam_spec, 'spec')

    # print summary of the results
    print("Teff: {:.2f} +/- {:.2f}".format(np.mean(flat_spec[:,0]), np.std(flat_spec[:,0])))
    print("logg: {:.2f} +/- {:.2f}".format(np.mean(flat_spec[:,1]), np.std(flat_spec[:,1])))
    print("[Fe/H]: {:.2f} +/- {:.2f}".format(np.mean(flat_spec[:,2]), np.std(flat_spec[:,2])))
    print("alpha: {:.2f} +/- {:.2f}".format(np.mean(flat_spec[:,3]), np.std(flat_spec[:,3])))

    # save the results to a csv file with the best fit parameters
    import pandas as pd

    results_df = pd.DataFrame()
    results_df['Teff'] = np.mean(flat_spec[:,0])
    results_df['logg'] = np.mean(flat_spec[:,1])
    results_df['feh'] = np.mean(flat_spec[:,2])
    results_df['alpha'] = np.mean(flat_spec[:,3])
    results_df['Teff_err'] = np.std(flat_spec[:,0])
    results_df['logg_err'] = np.std(flat_spec[:,1])
    results_df['feh_err'] = np.std(flat_spec[:,2])
    results_df['alpha_err'] = np.std(flat_spec[:,3])

    results_df.to_csv(os.path.join(plot_dir, 'results_spec.csv'), index=False)

