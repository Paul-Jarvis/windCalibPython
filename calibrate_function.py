# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 18:48:48 2022

@author: esnee
"""
import math

def calibrate(cam,P_pixel):
    
    # Camera details
    incl = math.radians(cam.incl)
    dist_diff = cam.dist2plane
    pixel_height = cam.pixel_height
    pixel_width = cam.pixel_width
    FOV_V = math.radians(cam.FOV_V)
    FOV_H = math.radians(cam.FOV_H)
    
    i = P_pixel[1]
    j = P_pixel[0]

    delta_theta_z = FOV_V / pixel_height
    delta_theta_h = FOV_H / pixel_width
    
    incl_h = math.radians(0)
    
    #convert all to radians
    
    # Calibrate for height
    
    diff_z = (dist_diff/2)*(math.tan(incl - (FOV_V/2) + (i-1)*delta_theta_z) + math.tan(incl - (FOV_V/2) + (i*delta_theta_z)))
    
    # Calibrate for length
    diff_h = (dist_diff/2)*(math.tan(incl_h - (FOV_H/2) + (j-1)*delta_theta_h) + math.tan(incl_h - (FOV_H/2) + (j*delta_theta_h)))

    
    return diff_z, diff_h