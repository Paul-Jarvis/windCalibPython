#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:15:46 2022

@author: pjarvis
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from calcInclination_function import calcInclination
from calibrate_function import calibrate
from extractWeather_function import extract_weather
from getAllUncertinity_function import getAllUncertinity

## Script to run windyPlume for example of Etna 2013

## %%%%%%%%%%%%%% User Input Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
## %%% Update pixels to claibrate %%%
## % - can be vectors for multiple points to calibrate
readPoints = 'n' #If 'y', read points from file (define dataFile, xCol and
                  #yCol). If 'n', define x_point, y_point.
# dataFile = '/scratch/paulj/sabancaya/processedData/20180808T1123/output/height.dat';
# xCol = 3 #Column of data file containing x-coordinate of pixels
# yCol = 8 #Column of data file containing y-coordinate of pixels
x_point = np.array([500, 600]) # x coordinate of pixel to claibrate in image frame
y_point = np.array([20, 30]) # y coordinate of pixel to claibrate in image frame

## %%% Set filename of the image frame %%%
imageFrame = 'exampleData/Etna2013/get_vent.jpg'
    
# %%% Set weather data %%%
# Vector of year, month, day, hour, minute, second of time at which to 
# calibrate
b = "2013-11-23 10:20"

#Names of netCDF files containing wind and geopotential data. If files are
#the same, give same name for both.
windFilename = 'exampleData/Etna2013/Nov_2013.nc'
geopotFilename = 'exampleData/Etna2013/Nov_2013.nc'

topPoint_Wind = 12000 # Height (in m a.s.l) which defines the upper limit 
                     #of height range to claulate the wind orientation over

## %%% Set camera properties %%%
class cam:
    pixel_width   = 704    # Width in pixels of image frame
    pixel_height  = 608    # Height in pixels of image frame
    FOV_H         = 18     # Horizontal field of view of the camera
    z_cam         = 35     # Height of the camera in m a.s.l
    oriCentreLine = 351.5  # Orientation of the camera to the centre of 
                           #the image frame
    dist2plane    = 26200  # Distance between the camera and the plane for 
                           #which points will be calibrated on
    ## %%% IF cam.incl in set to -100 i.e., unknown camera, inclination, update 
    ## %   cam.dist2ref, vent.z_ref and vent.centre_pixel_height %%%
    incl = -100.0 # Inclination of the camera in degrees (put as -100 if 
                # unknown)
    dist2ref             = 27300 # Distance between the camera loction and 
                                  #a known point in the image frame

min_FOVH = 16 # Minimum value of the horizontal field of view of the 
                   # camera (lower bound of the uncertinity)
max_FOVH = 20 # Maximum value of the horizontal field of view of the 
                   # camera (upper bound of the uncertinity)
    
class vent:
    z_ref               = 3300  # Height in m a.s.l of a known point in 
                                  #the image frame (same at the point 
                                  #defined on line 29)
    centre_pixel_height = 498   # Pixel value in the y direction of a 
                                  #known point in the image frame (same at 
                                  #the point defined on line 29)
    
## %%% Set vent properties
ventKnown = 'n'    #If 'n', needs to be determined
x_pixel_vent = 548 #Pixel coordinates of vent. Don't need to be entered 
                    #if ventKnown = 'n'
y_pixel_vent = 447
ventAlt = 5900     #Altitude of vent in m above sea level
ventLat = -15.75
ventLong = -71.75
    
## %%% Outfile
outFile = 'exampleData/Etna2013/windHeight.csv';
    
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## %% Calculate other camera properties

cam.FOV_V = cam.FOV_H *(cam.pixel_height/cam.pixel_width)     # Calculate vertical field of view of the camera 
    
if cam.incl == -100:
    incl = calcInclination(cam,vent) # Calculate inclination of the camera 
    cam.incl = incl
 
    
## %% Determine vent position
if ventKnown == 'n':
    A = plt.imread(imageFrame)   # Load in imageFrame

    plt.figure() # Display imageFrame
    plt.imshow(A)

    # Manually select the position of the vent from which the plume origintes from
    pts = plt.ginput(n=1,timeout=30, show_clicks=True)
    x_pixel_vent = pts[0][0]
    y_pixel_vent = pts[0][1]                              

# Put pixel position of the vent into vector where position (1) is the x
# position and (2) is the y position. y position is adjusted as 0 starts at
# the top 
P_vent = np.array([x_pixel_vent, cam.pixel_height-y_pixel_vent])
    
## %%% Calibrate vent position %%%
[diff_z, diff_h] = calibrate(cam, P_vent) #Turn pixel location to diff_z
                                           #(difference in calibrated point to
                                           #camera in the vertical) and diff_h
                                           #(difference in calibrated point to
                                           #camera in the horizontal)
    
vent_z = cam.z_cam + diff_z # Calulate absolute height (in m a.s.l) of the vent
vent_x = diff_h # Set vent_x to diff_h      

## %% Load weather data %%%

# Extract weather to a structure called netCDF
netCDF = extract_weather(windFilename, b, geopotFilename, ventLat, ventLong) 

## Get orientation  of the wind above the volcano %%%
wind_v_interp_f = interp1d(netCDF.z, netCDF.v) #Determine the wind velocity (v
                                              #component) of the range defined at
                                              #the base by the height of the vent
                                              #and at the top by the
                                              #topPoint_wind
nHeights = np.rint(topPoint_Wind - vent_z) + 1
nHeights = nHeights.astype(int)
z_array = np.linspace(vent_z, topPoint_Wind, nHeights)
wind_v = wind_v_interp_f(z_array)

wind_u_interp_f = interp1d(netCDF.z, netCDF.u) #Determine the wind velocity (u
                                              #component) of the range defined at
                                              #the base by the height of the vent
                                              #and at the top by the
                                              #topPoint_wind
wind_u = wind_u_interp_f(z_array)

ori = np.arctan2(wind_v, wind_u) #Calulate the wind oriention of the height
                                 #range of interest
ori = np.degrees(ori)

## Adjust angle of wind orientation to be with respect to the camera location %%%
Ori = ori
for i in range(nHeights):
    if ori[i] <= 90:
        Ori[i] = 90 - ori[i]
    else:
        Ori[i] = 360 - ori[i] + 90

## define Wind orientations to use in the calibration %%%
min_Ori = np.amin(Ori) # Determine min wind orientation of the range of the
                          # interest
max_Ori = np.amax(Ori) # Determine max wind orientation of the range of the
                          # interest
Ori = np.sum(Ori) / np.size(Ori) # Determine average wind orientation of the
                                 # range of the interest

## Read in points of interest
if readPoints == 'y':
    inData = np.genfromtxt(dataFile, delimiter = ',')
    x_point = inData[:, xCol - 1]
    y_point = cam.pixel_height - inData[:, yCol - 1]

##  Calibrate ponts of interest

# - height        == the height/s of the point/s of interest
# - lowerUncert_z == the minimum height/s of the point/s of interest
# - upperUncert_z == the maximum height/s of the point/s of interest
# - dist          == the difference in the horizontal poisiton/s with respect the the camera of the point/s of interest
# - lowerUncert_x == the minimum differne in the horizontal poisiton/s with respect the the camera of the point/s of interest
# - upperUncert_x == the maximum differne in the horizontal poisiton/s with respect the the camera of the point/s of interest
[upperUncert_x,lowerUncert_x,upperUncert_z,lowerUncert_z,dist,height] = getAllUncertinity(vent_x,vent_z,x_point,y_point,cam,Ori,P_vent,min_FOVH,max_FOVH,min_Ori,max_Ori);
     
## Need to output data    
np.savetxt(outFile, height);
