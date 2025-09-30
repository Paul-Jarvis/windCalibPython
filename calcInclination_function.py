# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 15:58:22 2022

@author: esnee
"""
import math

def calcInclination(cam,vent):
    #Camera details
    z_cam = cam.z_cam                                                          # height of camera in m's asl
    dist_diff = cam.dist2ref                                                   # distance between camera and reference point in m's asl
    pixel_height = cam.pixel_height
    FOV_V = cam.FOV_V
    
    z_ref = vent.z_ref
    
# Calculate Inclination
    phii = math.atan((z_ref-z_cam)/(dist_diff))
   # phii = math.degrees(phii)
    
    z_refPixel = pixel_height - vent.centre_pixel_height

    ## COMMENTED OUT LINE IN INCORREECT EXPRESSION, EQUATION 14 IN SNEE ET AL.
    #  (2023)
    #incl = phii +  math.atan((1 - ((2*z_refPixel)/pixel_height))*math.tan(math.radians(FOV_V/2)))
    incl = FOV_V * (0.5 - z_refPixel / pixel_height) + phii

    incl = math.degrees(incl)
    
    return incl
