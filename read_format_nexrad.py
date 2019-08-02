# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:17:40 2019
@author: Robinson Wallace

Purpose: To reorder nexrad data so that each scan has it's 0 azimuth degree 
         located in the 0th index for each radar elevation scan. PyArt does
         this behind the scenes, but we do this here for simplicity of analysis

"""

import os
import pyart
import numpy as np


def read_format_nexrad_lvl2(data_path):

    radexit = 0
    while True:
    
        # Make a reflectivity volume
        try:
            radar = pyart.io.read(data_path)
            radexit = 0
        except:
            # Sometimes there is an error in reading the radar file
            print('Something went wrong, retrying...')
            radexit =1
        
        if radexit == 1:
            print('Unable to load in NEXRAD file.')
            break
        
        rad_lat = radar.latitude['data'][0]
        rad_lon = radar.longitude['data'][0]
        radhght = radar.altitude['data'][0]
    
        sweep_start_inds = radar.sweep_start_ray_index['data']
        sweep_end_inds = radar.sweep_end_ray_index['data']
        n_azm = 360 #radar.nrays
        n_rng = radar.ngates
        n_elv = radar.nsweeps
        r_elv = np.zeros((n_elv))
        ref_vol = np.zeros((n_azm,n_rng,n_elv),dtype=np.float32)
        zdr_vol = np.zeros((n_azm,n_rng,n_elv),dtype=np.float32)
        rhv_vol = np.zeros((n_azm,n_rng,n_elv),dtype=np.float32)
        vel_vol = np.zeros((n_azm,n_rng,n_elv),dtype=np.float32)
        
        
        """
        The data comes as a long stream of data, not broken up into sweeps, so we
        have to do that. Some scans are 720 degrees while ohters are 360. We'll just
        drop the inbetween data for the 720 degree measurements. In the process, shift
        the azimuth scan so that north is at the first index
        """
            
        for e in range(0,n_elv):
            
            # For the sweeps of 720 azimuths
            if sweep_end_inds[e]-sweep_start_inds[e] == 719:
                f_azm = np.where(radar.azimuth['data'][sweep_start_inds[e]:sweep_end_inds[e]:2]<1)[0]
                # If the radar doesn't do a complete scan, the azimuths may not start at zero
                if len(f_azm) == 0:
                    continue
                ref = radar.fields['reflectivity']['data'][sweep_start_inds[e]:sweep_end_inds[e]:2,:]
                ref_sft = np.roll(ref,-f_azm[0],0)
                ref_vol[:,:,e] = ref_sft
                
                zdr = radar.fields['differential_reflectivity']['data'][sweep_start_inds[e]:sweep_end_inds[e]:2,:]
                zdr_sft = np.roll(zdr,-f_azm[0],0)
                zdr_vol[:,:,e] = zdr_sft
                
                rhv = radar.fields['cross_correlation_ratio']['data'][sweep_start_inds[e]:sweep_end_inds[e]:2,:]
                rhv_sft = np.roll(rhv,-f_azm[0],0)
                rhv_vol[:,:,e] = rhv_sft
                
                vel = radar.fields['velocity']['data'][sweep_start_inds[e]:sweep_end_inds[e]:2,:]
                vel_sft = np.roll(vel,-f_azm[0],0)
                vel_vol[:,:,e] = vel_sft
                
                r_elv[e] = radar.elevation['data'][sweep_start_inds[e]]
            
            if sweep_end_inds[e]-sweep_start_inds[e] == 359:
                f_azm = np.where(radar.azimuth['data'][sweep_start_inds[e]:sweep_end_inds[e]+1]<1)[0]
                if len(f_azm) == 0:
                    continue
                
                ref = radar.fields['reflectivity']['data'][sweep_start_inds[e]:sweep_end_inds[e]+1,:]
                ref_sft = np.roll(ref,-f_azm[0],0)
                ref_vol[:,:,e] = ref_sft
                
                zdr = radar.fields['differential_reflectivity']['data'][sweep_start_inds[e]:sweep_end_inds[e]+1,:]
                zdr_sft = np.roll(zdr,-f_azm[0],0)
                zdr_vol[:,:,e] = zdr_sft
                
                rhv = radar.fields['cross_correlation_ratio']['data'][sweep_start_inds[e]:sweep_end_inds[e]+1,:]
                rhv_sft = np.roll(rhv,-f_azm[0],0)
                rhv_vol[:,:,e] = rhv_sft
                
                vel = radar.fields['velocity']['data'][sweep_start_inds[e]:sweep_end_inds[e]+1,:]
                vel_sft = np.roll(vel,-f_azm[0],0)
                vel_vol[:,:,e] = vel_sft
                
                r_elv[e] = radar.elevation['data'][sweep_start_inds[e]]
                
        r_gate = radar.range['data']
    
        # Sometimes PyART doesn't decode the data correctly, so we can try it until it does
        # This has been coded to evacuate the loop if the data just truely sucks
        if len(r_elv) > 2:
            break
        else:
            print("Bad computation of elevation scans. Retrying...")
            radexit += 1
            if radexit == 2:
                break
    
    if radexit == 2:
        print("WARNING: Radar Data is erroneous (Bad elevation scans)")
        
    
    radar_time = radar.time['units'][14:]
    year = radar_time.split('-')[0]
    month = radar_time.split('-')[1]
    day = radar_time.split('-')[2][0:2]
    hour = radar_time.split('-')[2].split(':')[0][3:5]
    minute = radar_time.split('-')[2].split(':')[1]
    second = radar_time.split('-')[2].split(':')[2]
    
    if radexit >= 1: return None
    
    return {"ref":ref_vol,"zdr":zdr_vol,"rhv":rhv_vol,"velocity":vel_vol,
            "azim":np.arange(0,360),"gate":r_gate,"elv":r_elv,
            "year":year,"month":month,"day":day,"hour":hour,"minute":minute,
            "second":second,"radar_altitude":radhght,"radar_latitude":rad_lat,
            "radar_longitude":rad_lon}
    
if __name__ == '__main__':
    
    # Input a filepath to a nexrad data file
    vol = read_format_nexrad_lvl2('')
