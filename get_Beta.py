# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# Script to calibrate beta from using the plume shape and plume model

import numpy as np
import os
# import cv2
import imageio
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import math 
from calcInclination_function import calcInclination
from calibrate_function import calibrate
from extractWeather_function import extract_weather
from Calibrate_allTogether_function import Calibrate_allTogether
from getAllUncertinity_function import getAllUncertinity
# base_location = "H:\cameraCode"


x_select_ul = np.array([500, 600])
y_select_ul = np.array([20, 30])

b = "2013-11-23 10:20"    

w_filename = 'Nov_2013.nc'

# 
# define calibration constants
class cam:
    pixel_width = 704
    pixel_height = 608
    FOV_H = 18
    oriCentreLine = 180 + 171.5
    z_cam = 35
    dist2ref = 27300
    dist2plane = 26200
    FOV_V = FOV_H *(pixel_height/pixel_width)
    incl = 0
    
min_FOVH = 16
max_FOVH = 20


#cam = {'pixel_width':pixel_width, 'pixel_height':pixel_height, 'FOV_H':FOV_H, 'FOV_V':FOV_V, 'z_cam':35., 'oriCentreLine':oriCentreLine, 'dist2ref':27300., 'dist2plane':26200.}
#vent = {'z_ref':3300., 'centre_pixel_height':498.}

class vent:
    z_ref = 3300
    centre_pixel_height = 498

#Calculate iclination
# os.chdir('functions')
    
incl = calcInclination(cam,vent)

   
# cam['incl'] = incl
cam.incl = incl

#Determine values for the vent
# os.chdir('..')
 # 
im = plt.imread('get_vent.jpg')   
plt.figure()
plt.imshow(im)

pts = plt.ginput(n=1,timeout=30, show_clicks=True)
pts = pts[0]
P_vent =  np.array([pts[0], cam.pixel_height-pts[1]])
P_vent =  np.array([404.0296,106.6020])
plt.scatter(pts[0], pts[1])
plt.show()
#im = imageio.imread('get_vent.jpg')
# A = cv2.imread('get_vent.jpg');
    
    
# os.chdir('functions') 
[diff_z, diff_h] = calibrate(cam,P_vent)
      
vent_z = cam.z_cam + diff_z
vent_x = diff_h

# os.chdir('..')

# #  Extract weather
netCDF = extract_weather(w_filename,b)

z_array = np.linspace(vent_z,12000, num = 1000)


wind_v_interp_f = interp1d(netCDF.z,netCDF.v)#,z_array) #cam diff_z))
wind_v = wind_v_interp_f(z_array)

wind_u_interp_f = interp1d(netCDF.z,netCDF.u)#,z_array)#z_cam + diff_z))
wind_u = wind_u_interp_f(z_array)

for j in np.linspace(0,np.size(wind_u)-1,np.size(wind_u)):
    j = int(j)
    wind_v[j] = math.radians(wind_v[j])
    wind_u[j] = math.radians(wind_u[j])
    
ori = np.arctan2(wind_v,wind_u)

for j in np.linspace(0,np.size(wind_u)-1,np.size(wind_u)):
    j = int(j)
    ori[j] = math.degrees(ori[j])

Ori = []

for j in np.linspace(0,np.size(ori)-1,np.size(ori)):
    j = int(j)
    if ori[j] <= 90:
        Ori += [90 - ori[j]]
    else:
        Ori += [360-ori[j] + 90]
  

min_Ori = min(Ori);
max_Ori = max(Ori);
Ori = sum(Ori) /np.size(Ori)

# # # os.chdir('..')
y_select_ul =  cam.pixel_height- y_select_ul

if np.size(y_select_ul) == 1:#type(y_selecvt_ul) == float == int or type(y_select_ul) == int:
    y_select_ul = np.array([y_select_ul])
# # # #calibrate upper plume margin
# # # os.chdir(base_location + '/functions')
        
# # #Calibrate lower plume margin
[upperUncert,lowerUncert,upperUncert,lowerUncert, dist,height] = getAllUncertinity(vent_x,vent_z,x_select_ul,y_select_ul,cam,Ori,P_vent,min_FOVH,max_FOVH,min_Ori,max_Ori);

plt.figure()
plt.scatter(P_vent[0],P_vent[1])
plt.scatter(dist,height)
plt.show()


