# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 09:06:40 2019
@author: Robinson Wallace

Purpose: To compute the distance and bearing between two geographic points

Methods are taken from http://www.movable-type.co.uk/scripts/latlong.html
for the Distance and Bearing
"""

import numpy as np

def calcdistbear(lon1,lat1,lon2,lat2):
    
    """ Compute the Distance Betweeen Points """
    
    # Earth radius in meters
    Re = 6371000
    
    radn_lon1 = np.radians(lon1)
    radn_lat1 = np.radians(lat1)
    radn_lon2 = np.radians(lon2)
    radn_lat2 = np.radians(lat2)
    
    radn_del_lat = radn_lat2-radn_lat1
    radn_del_lon = radn_lon2-radn_lon1
    
    a = np.sin(radn_del_lat/2)**2 + (
               np.cos(radn_lat1) * np.cos(radn_lat2) * np.sin(radn_del_lon/2)**2)
    
    c = 2 * np.arctan2(a**0.5,(1-a)**0.5)
    
    d = Re*c
    
    
    """ Compute the azimuth between point """
    
    azim = np.arctan2( np.sin(radn_del_lon)*np.cos(radn_lat2) , 
                       (np.cos(radn_lat1)*np.sin(radn_lat2)) - 
                       (np.sin(radn_lat1) * np.cos(radn_lat2) * np.cos(radn_del_lon)))
    
    azim = (np.degrees(azim)+360) % 360
    
    # Round and reduce the number of decimals to print
    d_r = str(round(d, 2))
    azim_r = str(round(azim, 2))
    
    print('Distance: '+ str(round(d/1000.,2)) + f' kilometers / Bearing: {azim_r} degrees')
    
    return d,azim
    
    
if __name__ == '__main__':
    
    dist,azim = calcdistbear(-105,39,-106,40)
