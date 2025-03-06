#!/usr/bin/env python3
"""
Compute Alt-Az coordinates and airmass statistics for target stars during a given night.
"""

import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from astroplan import Observer
from astropy.table import Table
from tqdm import tqdm

def compute_star_airmass(observatory_lat, observatory_lon, observatory_height, date, stars, n_steps=1000):
    """
    Computes the minimum and median airmass for a list of stars over the night.

    Parameters:
    -----------
    observatory_lat : float
        Latitude of the observatory in degrees.
    observatory_lon : float
        Longitude of the observatory in degrees.
    observatory_height : float
        Height of the observatory in meters.
    date : str
        Date of observation in 'YYYY-MM-DD' format.
    stars : list of dicts
        List of target stars. Each star should be a dict with keys:
          - 'name': Name of the star.
          - 'ra': Right Ascension (e.g. '18h36m56.33635s').
          - 'dec': Declination (e.g. '+38d47m01.2802s').
    n_steps : int, optional
        Number of time steps to sample during the night (default is 1000).

    Returns:
    --------
    results : list of dicts
        A list containing a dictionary for each star with keys:
          - 'name': Star name.
          - 'min_airmass': Minimum airmass (np.nan if never above horizon).
          - 'median_airmass': Median airmass (np.nan if never above horizon).
          - 'sunset': Sunset time (ISO format).
          - 'sunrise': Sunrise time (ISO format).
    """
    # Setup the observatory location and observer
    location = EarthLocation(lat=observatory_lat*u.deg,
                             lon=observatory_lon*u.deg,
                             height=observatory_height*u.m)
    observer = Observer(location=location, timezone='UTC')

    # Define the time at midnight UTC on the given date
    time_midnight = Time(date + " 00:00:00")
    
    # Compute sunset and sunrise times for the night
    sunset = observer.sun_set_time(time_midnight, which='next')
    sunrise = observer.sun_rise_time(time_midnight, which='next')

    # Create a time grid spanning from sunset to sunrise
    night_times = sunset + (sunrise - sunset) * np.linspace(0, 1, n_steps)

    results = []
    for star in tqdm(stars):
        # Create a SkyCoord object for the target star (assumed to be in ICRS)
        star_coord = SkyCoord(star['ra']*u.deg, star['dec']*u.deg, frame='icrs')
        
        # Define the AltAz frame for the entire night at the observatory location
        altaz_frame = AltAz(obstime=night_times, location=location)
        
        # Transform the star's coordinates into the Alt-Az frame
        star_altaz = star_coord.transform_to(altaz_frame)
        
        # Compute airmass using the sec(z) attribute (only valid above horizon)
        airmass = star_altaz.secz
        
        # Consider only times when the star is above the horizon (altitude > 0)
        above_horizon = star_altaz.alt > 0*u.deg
        valid_airmass = airmass[above_horizon]
        
        if len(valid_airmass) > 0:
            min_airmass = np.min(valid_airmass)
            median_airmass = np.median(valid_airmass)
        else:
            min_airmass = np.nan
            median_airmass = np.nan
        
        results.append({
            'name': star['name'],
            'min_airmass': min_airmass,
            'median_airmass': median_airmass,
            'sunset': sunset.iso,
            'sunrise': sunrise.iso
        })
    
    return results

# Example usage
if __name__ == "__main__":
    # Define observatory parameters (LCO)
    observatory_lat = -29.015972222222      # degrees
    observatory_lon = -70.692083333333    # degrees
    observatory_height = 2380      # meters

    # Set the observation date
    date = "2025-03-13"


    # load catalog 
    data = Table.read('/Users/mncavieres/Documents/2024-1/Delve/filtered_stars.fits')

    # make a dictionary with source_id and coordinates
    stars = []
    for i in range(len(data)):
        stars.append({'name': data['source_id'][i], 'ra': data['ra'][i], 'dec': data['dec'][i]})


    # Compute the airmass statistics for each target star
    results = compute_star_airmass(observatory_lat, observatory_lon, observatory_height, date, stars)
    
    # save results to a csv using pandas
    import pandas as pd
    data_frame = pd.DataFrame(results)
    data_frame.to_csv('/Users/mncavieres/Documents/2024-1/Delve/Proposal/mike/airmass_results.csv', index=False)

    # merge results with the catalog
    input_data_df = data.to_pandas()
    results_df = pd.DataFrame(results)

    # merge the dataframes using the name column and source_id
    merged_df = pd.merge(input_data_df, results_df, left_on='source_id', right_on='name')
    merged_df.to_csv('/Users/mncavieres/Documents/2024-1/Delve/Proposal/mike/filtered_stars_with_airmass.csv', index=False)


    # Print the results
    for res in results:
        print(f"{res['name']}: Min Airmass = {res['min_airmass']:.2f}, Median Airmass = {res['median_airmass']:.2f}")
        print(f"   Sunset: {res['sunset']}, Sunrise: {res['sunrise']}")
