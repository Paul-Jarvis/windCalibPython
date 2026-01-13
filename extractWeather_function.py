# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 19:41:22 2022

@author: esnee
"""
import netCDF4
import numpy as np
from datetime import date
import sys
from scipy.interpolate import interp1d
from datetime import datetime as dt
import pandas as pd

def datenum(d):
    """
    MATLAB-like datenum for:
      - str (single timestamp)
      - pandas Timestamp
      - Python datetime (class aliased as dt)
      - pandas Series of datetimes (or strings parsable as datetimes)

    Returns:
      - float for scalar input
      - pandas Series of floats for Series input
    """
    # Normalize: strings → Timestamp; Series of strings → datetime64
    if isinstance(d, str):
        d = pd.to_datetime(d)  # scalar Timestamp
    elif isinstance(d, pd.Series):
        if not pd.api.types.is_datetime64_any_dtype(d):
            d = pd.to_datetime(d, errors='raise')  # vectorized parse

    # Scalar case: pandas Timestamp or Python datetime (aliased as dt)
    if isinstance(d, (pd.Timestamp, dt)):
        ordinal = d.toordinal()
        
        frac = (d - dt.fromordinal(ordinal)).total_seconds() / 86400.0
        return 366 + ordinal + frac

    # Series case
    if isinstance(d, pd.Series):
        # Ordinal per element: use date component
        ords = d.dt.date.apply(lambda date_obj: date_obj.toordinal())

        # Fractional day: seconds since midnight
        frac = (d - d.dt.floor('D')).dt.total_seconds() / 86400.0

        return 366 + ords + frac

    raise TypeError("datenum expects a str, pandas Timestamp, Python datetime, or Series of datetimes.")

class NetCDF:
    def __init__(self, z, u, v, v_a):
        self.z = z
        self.u = u
        self.v = v
        self.v_a = v_a

def extract_weather(windFile, b, geopotFile, vent_lat, vent_long):
    
# read in the netCDF file
    f_wind = netCDF4.Dataset(windFile)
    f_geopot = netCDF4.Dataset(geopotFile)
    
    z = f_geopot['z'][:]
    time = f_wind['time'][:]
    u = f_wind['u'][:]
    v = f_wind['v'][:]
    lat = f_geopot['latitude'][:]
    long = f_geopot['longitude'][:]
    
    z = np.squeeze(z)
    u = np.squeeze(u)
    v = np.squeeze(v)
    
    z = z/9.80665
    
    
    start_time = date.toordinal(date(1900,1,1))+366
    
    d = pd.to_datetime(b)
    end_time = datenum(d)

    time_want = (end_time - start_time) * 24

    if isinstance(time_want, pd.Series): #time_want is a sorted Pandas series
        if time_want.iloc[-1] > time[-1] or time_want[0] < time[0]:
            print('error with the selected times')
            sys.exit()
        else:
            if z.ndim == 4: #netCDF file contains data at multiple latitudes
                            #and/or longitudes

                # Determine indices of latitude and longitude needed in file
                row_lat = np.nonzero((lat - vent_lat) == 0)
                row_long = np.nonzero((long - vent_long) == 0)

                # Get reduced 2D arrays of weather data
                z = np.squeeze(z[:, :, row_lat, row_long])
                u = np.squeeze(u[:, :, row_lat, row_long])
                v = np.squeeze(v[:, :, row_lat, row_long])

            # Create interpolators
            z_interp_f = interp1d(time, np.transpose(z))
            u_interp_f = interp1d(time,np.transpose(u))#,time_want)
            v_interp_f = interp1d(time,np.transpose(v))#,time_want)

            # For each time, interpolate and find profiles of u(z) and v(z)
            new_z = z_interp_f(time_want)
            new_u = u_interp_f(time_want)
            new_v = v_interp_f(time_want)

            # Create output list and populate
            netCDF = []

            for i in range(len(time_want)):
                v_a = (new_u[:, i] ** 2 + new_v[:, i] ** 2) ** 0.5
                netCDF.append(NetCDF(new_z[:, i], new_u[:, i], new_v[:, i], v_a))
            
    else: #time_want is a single numeric value
        if time_want > time[-1] or time_want < time[0]:
            print('error with the selected time')
            sys.exit()

        ## select data for the time I want    
        if z.ndim == 4:
            ## In this case, the netCDF file contains data at multiple latitudes
            ## and/or longitudes
        
            row_t = np.nonzero((time - np.round(time_want)) == 0) # Select time in the 4D matrix
            row_lat = np.nonzero((lat - vent_lat) == 0) #Select latitude in the 4D matrix
            row_long = np.nonzero((long - vent_long) == 0) #Select longitude in the 4D matrix

            z = z[row_t, :, row_lat, row_long]
            u = u[row_t, :, row_lat, row_long]
            v = v[row_t, :, row_lat, row_long]

            new_z = np.squeeze(z)
            new_u = np.squeeze(u)
            new_v = np.squeeze(v)

        elif z.ndim == 2:
            ## In this case, the netCDF files contains data at a single latitude and
            ## longitude
            z_interp_f = interp1d(time,np.transpose(z))#,time_want)
            new_z = z_interp_f(time_want)

            u_interp_f = interp1d(time,np.transpose(u))#,time_want)
            new_u = u_interp_f(time_want)

            v_interp_f = interp1d(time,np.transpose(v))#,time_want)
            new_v = v_interp_f(time_want)        
        

        netCDF = NetCDF(new_z, new_u, new_v, (new_u ** 2 + new_v ** 2) ** 0.5)


    return netCDF
