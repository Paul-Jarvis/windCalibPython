#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 12:06:50 2022

@author: esnee
"""

import numpy as np
from calibrate_function import calibrate
from calibrateWind_function import calibrateWind

def Calibrate_allTogether(vent_x,vent_z,x_select,y_select,cam,Ori,P_vent):
    dist = np.zeros(np.size(x_select))
    height = np.zeros(np.size(x_select))

    for j in np.linspace(0,np.size(x_select)-1,np.size(x_select)):
        j = int(j)
        P_pixel = [x_select[j],y_select[j]]
        [diff_z, diff_x] = calibrate(cam,P_pixel)

        distanceFromVent_P1 = abs(diff_x - vent_x)

        [x,y,z,lambda_,w_tilde] = calibrateWind(Ori,cam,distanceFromVent_P1,P_vent,P_pixel)

        #print(z)
        #wait = input("Press Enter to continue.")
        
        if cam.oriCentreLine + 180 <360:
            omega_prime = cam.orCentrLine + 180
        elif cam.oriCentreLine + 180 >360:
            omega_prime = cam.oriCentreLine - 180
            
        if P_pixel[0] >= P_vent[0]:
            dist[j] = (((distanceFromVent_P1 + x)**2) + (y**2))**0.5
        elif P_pixel[0] < P_vent[0]:
            dist[j] = -(((distanceFromVent_P1 + x)**2) + (y**2))**0.5
            
        if w_tilde > 180 and w_tilde < 360:
            dist[j] = dist[j] * -1
            
        height[j] = cam.z_cam +diff_z + z
            
    return dist, height, lambda_
        
       



