


#/////////////////////
#  open_buoy_data ///
#///////////////////
#---------------------------------------------------------------------
# Opens and returns time-cropped SIDEx buoy coordinates
#---------------------------------------------------------------------
# DEPENDENCIES
import numpy as np
import pandas as pd
import datetime
#--------------------------------------------------------------------

def open_buoy_data(buoy_file, 
                   start_date = datetime.datetime(year=2021,month=1,day=1,hour=0),
                   end_date = datetime.datetime(year=2022,month=1,day=1,hour=0)):
    
    """Open and return time-cropped SIDEx buoy coordinates.
    Updated version for new version of buoy track data produced by Jenny/Angela Jan 2022.
    Data should have columns in order: ['latitude', 'longitude', 'datetime'].

INPUT: 
- buoy_file: buoy file from  SIDEx buoy track data produced by Jenny/Angela Jan 2022.
            (e.g. path+'OSU-IT-53_300534061988500.csv')
- start_date: datetime object corresponding to desired initial date of buoy data to import
- end_date: datetime object corresponding to desired final date of buoy data to import
    * note: if wanting to import all data from buoy file, simply provide start_date/end_date much outside date range of buoy data. default start/end dates are outside SIDEx period, so by default imports all buoy data.

OUTPUT:
- lon_track: array of buoy longitudes between start/end date
- lat_track: array of buoy latitudes between start/end date
- time_track: array of times corresponding to coordinates
- df_buoy: pandas data frame of opened and time-cropped data

DEPENDENCIES:
import numpy as np
import pandas as pd
import datetime

Latest recorded update:
07-20-2023
    """
    
    # read in data and convert datetime column to datetime objects
    # this method is specific to the way the SIDEx buoy data are stored
    df_buoy = pd.read_csv(buoy_file, parse_dates=['datetime'],
                          date_parser = lambda x: datetime.datetime.strptime(x, "%Y-%m-%d %H:%M:%S"))
    
    # make sure dates are in datetime format
    df_buoy['datetime'] = pd.to_datetime(df_buoy['datetime'])  
    
    # select only dates that fall between/including start and end dates
    mask = (df_buoy['datetime'] >= start_date) & (df_buoy['datetime'] <= end_date)
    df_buoy = df_buoy.loc[mask]

    # make sure at least some dates fall within range, else return empty data
    if len(df_buoy) > 0:
        lon_track = df_buoy.longitude.values
        lat_track = df_buoy.latitude.values
        time_track = pd.to_datetime(df_buoy.datetime.values)

    # if no dates within range, return empty list
    else:
        lon_track = []
        lat_track = []
        time_track = []
        
    return lon_track, lat_track, time_track, df_buoy



from metpy.units import units
import numpy as np
import pandas as pd
import datetime

# Import the geodesic module from geopy library 
from geopy.distance import geodesic
from pyproj import Geod
g = Geod(ellps='WGS84')




def calc_velocity(lon_track = [], lat_track = [], time_track = [], step = 1, 
                  vflag = 100*units('cm')/units('s')):

    
#     # step size between buoy data, usually step = 1 is 10 minutes
#     # but some steps larger
#     #====================
#     step = 2
#     #====================

# vflag: velocity flag, mask data with velocity components greater than this

    # times halfway between considered time steps
    time = np.array([])

    # grab every nth step
    lons = lon_track[::step]
    lats = lat_track[::step]
    times = time_track[::step]


    distance = np.array([])
    azimuth = np.array([])
    u = np.array([])
    v = np.array([])
    dx = np.array([])
    dy = np.array([])
    sp = np.array([])
    dt = np.array([])

    for ii in range(len(lons)-1):

        # start location
        loc1=(lats[ii], lons[ii])
        # end location
        loc2=(lats[ii+1], lons[ii+1])

        # calculate time length
        DT = (times[ii+1] - times[ii]).total_seconds() * units('second')
        dt = np.append(dt, DT)

        # save time between steps
        mid_time = times[ii] + datetime.timedelta(seconds = 0.5*DT.magnitude)
        time = np.append(time, mid_time)

        # compute forward and back azimuths, plus distance
        az12,az21,dist = g.inv(loc1[1],loc1[0],loc2[1],loc2[0])
        DI = dist*units('meter').to('cm')
        distance = np.append(distance, DI)
        azimuth  = np.append(azimuth, az12)

        # angle from east
        beta = 90 * units('degree') - az12 * units('degree')
        if beta <= -180*units('degree'):
            beta += 360*units('degree')

        # calculate zonal, meridional displacements fmor azimuth
        DX = (dist*units('meter') * np.cos(beta.to('radian'))).to('cm')
        DY = (dist*units('meter') * np.sin(beta.to('radian'))).to('cm')
        
        dx = np.append(dx, DX)
        dy = np.append(dy, DY)
        u = np.append(u, DX/DT)
        v = np.append(v, DY/DT)
        sp = np.append(sp, DI/DT)
        
        
    # flag unlikely velocity values
    u[np.abs(sp) > vflag] = np.nan
    v[np.abs(sp) > vflag] = np.nan
    sp[np.abs(sp) > vflag] = np.nan

    return u, v, sp, time, dx, dy, distance