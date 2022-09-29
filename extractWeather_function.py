# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 19:41:22 2022

@author: esnee
"""
import netCDF4
import numpy as np
from datetime import date
from scipy.interpolate import interp1d
from datetime import datetime as dt

def datenum(d):
    return 366+ d.toordinal() + (d - dt.fromordinal(d.toordinal())).total_seconds()/(24*60*60)



def extract_weather(filename,b):
    
# read in the netCDF file
    f = netCDF4.Dataset(filename)
    
    z = f.variables['z']
    time = f.variables['time']
    pres = f.variables['level']
    temp = f.variables['t']
    u = f.variables['u']
    v = f.variables['v']
    r_h = f.variables['r']
    
    temp = np.squeeze(temp)
    z = np.squeeze(z)
    u = np.squeeze(u)
    v = np.squeeze(v)
    r_h = np.squeeze(r_h)
    
    z = z/9.80665
    
    
    start_time = date.toordinal(date(1900,1,1))+366
    
    d = dt.strptime(b,'%Y-%m-%d %H:%M')
    end_time =  datenum(d)
    
    time_want = (end_time - start_time)*24
    time = np.squeeze(time)
    
    temp_interp_f = interp1d(time,np.transpose(temp))#,time_want)
    new_temp = temp_interp_f(time_want)
    
    z_interp_f = interp1d(time,np.transpose(z))#,time_want)
    new_z = z_interp_f(time_want)
    
    u_interp_f = interp1d(time,np.transpose(u))#,time_want)
    new_u = u_interp_f(time_want)
    
    v_interp_f = interp1d(time,np.transpose(v))#,time_want)
    new_v = v_interp_f(time_want)
    
    class netCDF:
        temp = new_temp
        z = new_z
        u = new_u
        v = new_v
        v_a = (u**2 + v**2)**0.5


    return netCDF




# a = extract_weather(filename,b)