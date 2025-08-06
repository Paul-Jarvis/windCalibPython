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

def datenum(d):
    return 366+ d.toordinal() + (d - dt.fromordinal(d.toordinal())).total_seconds()/(24*60*60)



def extract_weather(windFile, b, geopotFile, vent_lat, vent_long):
    
# read in the netCDF file
    f_wind = netCDF4.Dataset(windFile)
    f_geopot = netCDF4.Dataset(geopotFile)
    
    z = f_geopot['z'][:]
    time = f_wind['time'][:]
    #pres = f.variables['level']
    #temp = f.variables['t']
    u = f_wind['u'][:]
    v = f_wind['v'][:]
    #r_h = f.variables['r']
    lat = f_geopot['latitude'][:]
    long = f_geopot['longitude'][:]
    
    #temp = np.squeeze(temp)
    z = np.squeeze(z)
    u = np.squeeze(u)
    v = np.squeeze(v)
    #r_h = np.squeeze(r_h)
    
    z = z/9.80665
    
    
    start_time = date.toordinal(date(1900,1,1))+366
    
    d = dt.strptime(b,'%Y-%m-%d %H:%M:%S')
    end_time =  datenum(d)
    
    time_want = (end_time - start_time)*24
    #time = np.squeeze(time)
    #temp_interp_f = interp1d(time,np.transpose(temp))#,time_want)
    #new_temp = temp_interp_f(time_want)

    if time_want > time[-1] or time_want < time[1]:
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
        

    class netCDF:
        #temp = new_temp
        z = new_z
        u = new_u
        v = new_v
        v_a = (u**2 + v**2)**0.5


    return netCDF




# a = extract_weather(filename,b)
