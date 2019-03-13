# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 11:21:33 2019
@author: Robinson Wallace - robinsonwwallace@gmail.com

Purpose: This routine parses out the data obtained from:
         http://apollo.lsc.vsc.edu/classes/remote/lecture_notes/radar/88d/88D_locations.html
         to get the specified radar's latitude and longitude. The source data
         gives the location in deg,min,sec format, but ultimate it's nice to
         have that in decminal format. This routine does that as well.
"""


import numpy as np
import csv

def radar_lonlat(radar_id,radarlocations_filepath):
    
    with open(radarlocations_filepath, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if row[1] == radar_id:
                
                latitude_degminsec  = row[3].split(' / ')[0]
                latitude_deg = np.float(latitude_degminsec[0:2])
                latitude_min = np.float(latitude_degminsec[2:4])
                latitude_sec = np.float(latitude_degminsec[4:])
                
                longitude_degminsec = row[3].split(' / ')[1]
                longitude_deg = np.float(longitude_degminsec[0:3])
                longitude_min = np.float(longitude_degminsec[3:5])
                longitude_sec = np.float(longitude_degminsec[5:])
                
                latitude_decimal = (latitude_deg) + (latitude_min/60.) + (latitude_sec/3600.)
                longitude_decimal = -1 * ((longitude_deg) + (longitude_min/60.) + (longitude_sec/3600.))
                
                print('Radar Latitude: '+str(round(latitude_decimal,2)), end=' / ')
                print('Radar Longitude: '+str(round(longitude_decimal,2)))
    
    return longitude_decimal, latitude_decimal
                
                
if __name__ == '__main__':
    
    radarlocations_filepath = 'C:/Users/robin/OneDrive/Graduate School/ATOC 7500 - Snow Observations/Python Programs/SnowGuage Radar Comparison/radarlocations.csv'
    
    location = radar_lonlat('KFTG',radarlocations_filepath)
