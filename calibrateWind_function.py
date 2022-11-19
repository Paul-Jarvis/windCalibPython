#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 07:22:57 2022

@author: esnee
"""
import math


def calibrateWind(Ori,cam,distanceFromVent,P_vent,P_pixel):
    incl = cam.incl
    FOV_V = cam.FOV_V
    FOV_H = cam.FOV_H
    pixel_height = cam.pixel_height
    pixel_width = cam.pixel_width
    oriCentreLine = cam.oriCentreLine
    
    if Ori <180:
        w_prime = Ori
    elif Ori >= 180:
        w_prime = Ori - 180
    else:
        print('Error with w_prime setting')
        
    if 0 < oriCentreLine and oriCentreLine <90:
        I = oriCentreLine + 90
    elif 90 < oriCentreLine and oriCentreLine <=270:
        I = oriCentreLine - 90
    elif 270 < oriCentreLine and oriCentreLine < 360:
        I = oriCentreLine -270
    else:
        print('Error with setting wind orientation')
        
        
    if abs(I - w_prime) < 90:
        lambda_ = abs(I - w_prime)
    elif abs(I - w_prime) > 90:
        lambda_ = 180 - abs(I - w_prime)
    else:
        print('Error with setting lambda')
        
    if (Ori - oriCentreLine) > 0:
        w_tilde = Ori - oriCentreLine
    elif Ori - oriCentreLine < 0:
        w_tilde = Ori - oriCentreLine + 360
    else:
        print('Error with setting w_tilde')


    alpha = P_pixel[0] * (FOV_H/pixel_width)     
    b = distanceFromVent
    chi = (P_pixel[1] * FOV_V)/ pixel_height
    di = incl - (FOV_V/2) + chi

    if (90 < w_tilde) and (w_tilde < 180) or (270 <= w_tilde) and (w_tilde < 360):
        cos_h = math.cos(math.radians(alpha - (FOV_H/2) - lambda_))
        sin_lambda = math.sin(math.radians(lambda_))
        h = (b * sin_lambda) / cos_h
        
        if P_pixel[0] >= P_vent[0]:
            x = -h * math.sin(math.radians((FOV_H/2)-alpha))
            y = h * math.cos(math.radians((FOV_H/2)-alpha))
            z = y * math.tan(math.radians(di))
        elif P_pixel[0] < P_vent[0]:
            x = h * math.sin(math.radians((FOV_H/2)-alpha))
            y = -h * math.cos(math.radians((FOV_H/2)-alpha))
            z = y * math.tan(math.radians(di))
            
    elif (0 < w_tilde) and (w_tilde < 90) or (180 <= w_tilde) and (w_tilde < 270):
        cos_h = math.cos(math.radians(alpha - (FOV_H/2) + lambda_))
        sin_lambda = math.sin(math.radians(lambda_))
        h = (b * sin_lambda) / cos_h
        
        if P_pixel[0] >= P_vent[0]:
            x = -h * math.sin(math.radians((FOV_H/2)-alpha))
            y = h * math.cos(math.radians((FOV_H/2)-alpha))
            z = y * math.tan(math.radians(di))
        elif P_pixel[0] < P_vent[0]:
            x = h * math.sin(math.radians((FOV_H/2)-alpha))
            y = -h * math.cos(math.radians((FOV_H/2)-alpha))
            z = y * math.tan(math.radians(di))
            
            
    return x, y, z, lambda_, w_tilde
        
        
        
        
        
        
