# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 18:44:51 2019
@author: Robinson Wallace - walacer@colorado.edu

Purpose: This routine creates a timeseries of reflectivity and differential
         reflectivity at the chosen location (given by azimuth and range
         from the specified radar). Data are output into a CSV file.
    
Required Modules: Most modules are co-installed with Anaconda. The user will 
                  need: nexradaws (https://pypi.org/project/nexradaws/) and 
                  read_format_nexrad_lvl2 (custom function: Robinson Wallace). 
                  
                  Both read_format_nexrad_lvl2 and calcdistbear are custom 
                  modules that should be distribtuted alongside this routine. 
                  If the user needs either module, contact Robinson for a copy.
                  
"""
import os
import csv
import pytz
import datetime
import tempfile
import numpy as np
from calcdistbear import calcdistbear
from read_format_nexrad import read_format_nexrad_lvl2
from radar_lonlat import radar_lonlat
import nexradaws

# Finds the index number of the value in the array closest to the desired value
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

# Converts data from log-space to linear-space and computes the averages 
# across input values before finally returning data back to log-space
def log_avg(values):
    values = np.array(values)
    no_logs = 10.0**(values/10.0)
    no_logs_mean = np.mean(no_logs)
    log_avg = 10*np.log10(no_logs_mean)
    return log_avg

""" ####################################################################### """
""" ############################# USER INPUTS ############################# """

# Designate the radar you want to retieve data for
radar_id = 'KFTG'

# Gauge longtidue and Latitude. Can be obtained by using Google Maps or 
# information provided by the station data
glon,glat = -105.242117, 40.010494

# Set the beginning and end times for the period of analysis in YYYY-MM-DD-HH:MI format
radar_begin_time = datetime.datetime.strptime('2019-03-03-20:00', '%Y-%m-%d-%H:%M')
radar_end_time   = datetime.datetime.strptime('2019-03-03-20:55', '%Y-%m-%d-%H:%M')

# The elevation scan to assess data at. Note: this is an integer value and
# not the actual degree value of the elevation scan. So, 0 is the lowest
# scan, 1 is the second lowest, etc. For gauge comparisons to the radar
# data, the lowest scan (0) is recommended.
elv = 0

# The failepath to the csv that contains the radar locations. These data
# can be obtained from http://apollo.lsc.vsc.edu/classes/remote/lecture_notes/radar/88d/88D_locations.html
# and need to be copied into a csv file.
radarlocations_filepath = '/path/to/radarlocations.csv'

# Where to save the generated CSV file
csv_filepath = '/path/to/csv_out_file.csv'

""" ####################################################################### """
""" ################## GET THE NEXRAD RADAR DATA FROM AWS ################# """

if radar_begin_time > radar_end_time:
    raise ValueError("Timespan for search invalid. End time needs to be after the start time")

# You can change the time zone if it's more conveient for you.
timezone = pytz.timezone('UTC')

# Set the begin and end times to a format the nexrad data call funtion can use
start = timezone.localize(radar_begin_time)
end = timezone.localize(radar_end_time)

templocation = tempfile.mkdtemp()
conn = nexradaws.NexradAwsInterface()

# This will scan the AWS servers for data between the specified time frame
scans = conn.get_avail_scans_in_range(start, end, radar_id)
if len(scans) == 0:
    raise ValueError("There are {} scans available between: \n{} and {}\n### Try expanding your search timeframe ###".format(len(scans), start, end))
else:
    print("There are {} scans available between {} and {}\n".format(len(scans), start, end))

# This will download the data to a temporary folder. On Windows this is in:
# C:/Users/<user>/AppData/Local/Temp/<randomly generated temp folder name>
results = conn.download(scans, templocation)
print(' ')

# Initialize empty lists 
ref_vector     = ['Reflectivity (dBZ)']
ref_avg_vector = ['Average Reflectivity (dBZ)']
zdr_vector     = ['Differential Reflectivity (dB)']
zdr_avg_vector = ['Average Differential Reflectivity (dB)']
rhv_vector     = ['Cross Correlation Coefficient']
vel_vector     = ['Velocity (m/s)']

yer_vec = ['Year']
mon_vec = ['Month']
day_vec = ['Day']
hor_vec = ['Hour']
min_vec = ['Minute']
sec_vec = ['Second']

# Get the radar latitude and longitude based on the input radar_id
rlon,rlat=radar_lonlat(radar_id,radarlocations_filepath)
print('Gauge Latitude: '+ str(round(glat,2)) +' / Gauge Longitude: '+str(round(glon,2)))

# Get the distance and azimuth between the radar and the gauge location
point_dist, point_azim = calcdistbear(rlon,rlat,glon,glat)
print(' ')

for scan in results.iter_success():
    
    print('Reading in: '+scan.filename)
    
    vol = read_format_nexrad_lvl2(scan.filepath)
    
    if vol is None:
        print('Unable to decode current volume. Skipping volume')
        continue
    
    cur_time = datetime.datetime.strptime(vol['year']+vol['month']+vol['day']+
                                          vol['hour']+vol['minute'],'%Y%m%d%H%M')
    
    # Ignore any data downloaded that isn't in the indicated timeframe. This
    # is most useful if the code is written to include a manually downloaded
    # large dataset. 
    if cur_time < radar_begin_time: continue
    if cur_time > radar_end_time: break

    # Find the nearest azimuth and gate index (since these need to be integers)
    azim = find_nearest(vol['azim'],point_azim)
    gate = find_nearest(vol['gate'],point_dist)
    
    # If the user puts in an elevation scan that exceeds the number of 
    # elevations, sidestep the resultant index error by reducing the 
    # scan height index to the second lowest (so the average can take
    # the gate above as well).
    if elv >= len(vol['elv']): 
        print('Chosen Elevation Scan was too large. Selecting second lowest scan.')
        print(' ')
        elv = len(vol['elv'])-2
    
    # Fill in the exact reflectivity and differential reflectivity at the gate 
    # closest to the input location and at the lowest elevation scan
    ref_vector.append(vol['ref'][azim,gate,elv])
    zdr_vector.append(vol['zdr'][azim,gate,elv])
    rhv_vector.append(vol['rhv'][azim,gate,elv])
    vel_vector.append(vol['velocity'][azim,gate,elv])
    
    # We can compute an average of reflectivity using the adjacent gates
    # This may provide more realistic values for the precip that reaches
    # the ground. Note: this does not account for temporal differences
    # in the time radar data and the ground data are respectively observed.
    ref_avg_vector.append(log_avg([vol['ref'][azim-1,gate,elv],
                                  vol['ref'][azim+1,gate,elv],
                                  vol['ref'][azim,gate-1,elv],
                                  vol['ref'][azim,gate+1,elv],
                                  vol['ref'][azim,gate,elv+1]]))
    
    
    """ There is something funky about averaging the ZDR values. """
    zdr_avg_vector.append(log_avg([vol['zdr'][azim-1,gate,elv],
                                  vol['zdr'][azim+1,gate,elv],
                                  vol['zdr'][azim,gate-1,elv],
                                  vol['zdr'][azim,gate+1,elv],
                                  vol['zdr'][azim,gate,elv+1]]))
    
    yer_vec.append(vol['year'])
    mon_vec.append(vol['month'])
    day_vec.append(vol['day'])
    hor_vec.append(vol['hour'])
    min_vec.append(vol['minute'])
    sec_vec.append(vol['second'][0:2])


# Write out the computed data to a csv file
with open(csv_filepath, mode='w', newline='') as targ:
    
    writer = csv.writer(targ)
    writer.writerows([yer_vec,
                      mon_vec,
                      day_vec,
                      hor_vec,
                      min_vec,
                      sec_vec,  
                      ref_vector,
                      ref_avg_vector,
                      zdr_vector,
                      zdr_avg_vector,
                      rhv_vector,
                      vel_vector])
    
