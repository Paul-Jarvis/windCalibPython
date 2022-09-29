#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:27:46 2022

@author: esnee
"""

from Calibrate_allTogether_function import Calibrate_allTogether
import numpy as np
import copy

def getAllUncertinity(vent_x,vent_z,x_select,y_select,cam, Ori, P_vent, min_FOVH,max_FOVH, min_Ori,max_Ori):
    
    [dist, height, lambda_] = Calibrate_allTogether(vent_x,vent_z,x_select,y_select,cam,Ori,P_vent)
    
    
    x_select_upper = np.zeros(np.size(x_select))
    x_select_lower = np.zeros(np.size(x_select))
    y_select_upper = np.zeros(np.size(x_select))
    y_select_lower = np.zeros(np.size(x_select))
    
    x_cal_camRes_min = np.zeros(np.size(x_select))
    x_cal_camRes_max = np.zeros(np.size(x_select))
    y_cal_camRes_min = np.zeros(np.size(x_select))
    y_cal_camRes_max = np.zeros(np.size(x_select))
    
    #Uncertinity from pixel 
    for l in np.linspace(0,np.size(x_select)-1,np.size(x_select)):
        l = int(l)
        x_select_upper[l] = x_select[l]+1
        x_select_lower[l] = x_select[l]-1
        y_select_upper[l] = y_select[l]+1
        y_select_lower[l] = y_select[l]-1
        
        #horizontal
    [dist_a, height_a, lambda_a] = Calibrate_allTogether(vent_x,vent_z,x_select_upper,y_select,cam,Ori,P_vent)
    
    [dist_b, height_b, lambda_a] = Calibrate_allTogether(vent_x,vent_z,x_select_lower,y_select,cam,Ori,P_vent)
    
    for l in np.linspace(0,np.size(x_select)-1,np.size(x_select)):
        l = int(l)
        o = abs((dist-dist_a)/2)
        p = abs((dist - dist_b)/2)
        x_cal_camRes_min[l] = min(o[l],p[l])
        x_cal_camRes_max[l] = max(o[l],p[l])
        
    
    
    #     #vertical
    [dist_c, height_c, lambda_c] = Calibrate_allTogether(vent_x,vent_z,x_select,y_select_upper,cam,Ori,P_vent)
    
    [dist_d, height_d, lambda_d] = Calibrate_allTogether(vent_x,vent_z,x_select,y_select_lower,cam,Ori,P_vent)
    
    
    for l in np.linspace(0,np.size(x_select)-1,np.size(x_select)):
        l = int(l)
        o = abs((dist-dist_c)/2)
        p = abs((dist - dist_d)/2)
        y_cal_camRes_min[l] = min(o[l],p[l])
        y_cal_camRes_max[l] = max(o[l],p[l])
    
    
    
    # #Uncertinity from FOV
    #     #Min FOV MAY NOT WORK cam_1 = CAM 
    cam_1 = copy.copy(cam)
    cam_1.FOV_H = min_FOVH
    cam_1.FOV_V = cam_1.FOV_H * (cam_1.pixel_height/cam_1.pixel_width)
    
    [dist_1, height_1, lambda_1] = Calibrate_allTogether(vent_x,vent_z,x_select,y_select,cam_1,Ori,P_vent)
    
    #     #Max FOV
    cam_2 = copy.copy(cam)
    cam_2.FOV_H = max_FOVH
    cam_2.FOV_V = cam_2.FOV_H * (cam_2.pixel_height/cam_2.pixel_width)
    
    [dist_2, height_2, lambda_2] = Calibrate_allTogether(vent_x,vent_z,x_select,y_select,cam_2,Ori,P_vent)
    
    #     #Determine errors
    x_cal_FOV_min = []
    x_cal_FOV_max = []
    
    z_cal_FOV_min = []
    z_cal_FOV_max = []
    
    for k in np.linspace(0,np.size(x_select)-1,np.size(x_select)):
          k = int(k)
          x_cal_FOV_min += [min(dist_1[k], dist_2[k])]
          x_cal_FOV_max += [max(dist_1[k], dist_2[k])]
         
          z_cal_FOV_min += [min(height_1[k], height_2[k])]
          z_cal_FOV_max += [max(height_1[k], height_2[k])]
         
      #Uncertinity from wind
          #min wind ori
    [dist_3, height_3, lambda_3] = Calibrate_allTogether(vent_x,vent_z,x_select,y_select,cam,min_Ori,P_vent)
        
        #max wind ori
    [dist_4, height_4, lambda_4] = Calibrate_allTogether(vent_x,vent_z,x_select,y_select,cam,max_Ori,P_vent)
    
        #determine errors
    x_cal_wind_min = []
    x_cal_wind_max = []
    
    z_cal_wind_min = []
    z_cal_wind_max = []
    
    for k in np.linspace(0,np.size(x_select)-1,np.size(x_select)):
          k = int(k)
          x_cal_wind_min += [min(dist_3[k], dist_4[k])]
          x_cal_wind_max += [max(dist_3[k], dist_4[k])]
         
          z_cal_wind_min += [min(height_3[k], height_4[k])]
          z_cal_wind_max += [max(height_3[k], height_4[k])]
          
          
    upperUncert_x = x_cal_FOV_max
    upperUncert_z = z_cal_FOV_max
    
    lowerUncert_x = x_cal_FOV_min
    lowerUncert_z = z_cal_FOV_min
    
    return upperUncert_x,lowerUncert_x,upperUncert_z,lowerUncert_z, dist, height
         
         
         
         
         
          
         
         