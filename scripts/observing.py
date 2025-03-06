#!/usr/bin/env python3
"""
Script to compute Alt–Az coordinates for a list of stars from a given observatory and date,
and to determine the minimum and median airmass for each star during the night.
"""

import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from astroplan import Observer

# ---------------------------
# Input observatory parameters
# ---------------------------
# Example: Mauna Kea Observatory
observatory_lat = 19.8207      # degrees
observatory_lon = -155.4681    # degrees
observatory_height = 4205      # in meters

# Create an EarthLocation object for the observatory
location = EarthLocation(lat=observatory_lat*u.deg,
                         lon=observatory_lon*u.deg,
                         height=observatory_height*u.m)

# Create an Observer object (astroplan uses this for computing sunset/sunrise)
observer = Observer(location=location, timezone='UTC')

# ---------------------------
# Define the observing date
# ---------------------------
# Example: Night of March 13, 2025.
date = '2025-03-13'
# Using midnight UTC for the date (adjust if needed for your local time)
time_midnight = Time(date + " 00:00:00")

# Compute sunset (after midnight) and sunrise times for that night
sunset = observer.sun_set_time(time_midnight, which='next')
sunrise = observer.sun_rise_time(time_midnight, which='next')

print(f"Sunset: {sunset.iso}")
print(f"Sunrise: {sunrise.iso}")

# ---------------------------
# Create a time grid for the night
# ---------------------------
n_steps = 1000
night_times = sunset + (sunrise - sunset) * np.linspace(0, 1, n_steps)

# ---------------------------
# List of stars with RA and Dec (ICRS)
# ---------------------------
stars = [
    {'name': 'Vega', 'ra': '18h36m56.33635s', 'dec': '+38d47m01.2802s'},
    {'name': 'Betelgeuse', 'ra': '05h55m10.30536s', 'dec': '+07d24m25.4304s'},
    # You can add more stars here
]

# ---------------------------
# Loop over each star to compute Alt-Az and airmass statistics
# ---------------------------
print("\nStar observations:")
for star in stars:
    # Create a SkyCoord object for the star
    star_coord = SkyCoord(star['ra'], star['dec'], frame='icrs')
    
    # Create an AltAz frame for the entire night
    altaz_frame = AltAz(obstime=night_times, location=location)
    
    # Transform the star's ICRS coordinates into Alt-Az coordinates
    star_altaz = star_coord.transform_to(altaz_frame)
    
    # Compute airmass using the sec(z) attribute.
    # Note: airmass is only meaningful when the star is above the horizon.
    airmass = star_altaz.secz
    
    # Create a mask for times when the star is above the horizon (altitude > 0 degrees)
    above_horizon = star_altaz.alt > 0*u.deg
    valid_airmass = airmass[above_horizon]
    
    if len(valid_airmass) > 0:
        min_airmass = np.min(valid_airmass)
        median_airmass = np.median(valid_airmass)
    else:
        min_airmass = np.nan
        median_airmass = np.nan
    
    print(f"{star['name']}: Min Airmass = {min_airmass:.2f}, Median Airmass = {median_airmass:.2f}")

    # Optionally, if you want to inspect the Alt-Az track, you could plot or print the coordinates
