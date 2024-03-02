# import the commonly used packages to save time see world

# system
import sys
import os
import glob

# math
import numpy as np
import numpy.ma as ma

# file management
import xarray as xr
import pandas as pd

# time
from datetime import datetime, timedelta
import calendar

# plotting
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.pyplot as plt
import cmocean

# metpy
import metpy.calc
from metpy.units import units

# geometry
import shapely
from shapely import wkt
from shapely.geometry import LineString
# specify a named ellipsoid
from pyproj import Geod
geod = Geod(ellps="WGS84")

# ignore warnings related to geographic plotting
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 

# homemade functions
sys.path.append('./scripts/')

from LIB_geo_func import *
from LIB_geo_plot import *
from LIB_plotting import *
from LIB_PPdrift_NSIDC0116 import *
from LIB_SIC_MODISAMSR2 import *
from LIB_access_ERA5 import *
from LIB_SIDExbuoy import calculate_velocity, open_buoy_data

# from LIB_geo_func import make_polygon, add_points_to_segment, within_polygon_indices
# from LIB_geo_plot import (add_land, add_coast, add_grid, add_date, fix_cartopy_vectors)
# from LIB_plotting import add_colorbar

# from LIB_PPdrift_NSIDC0116 import grab_ice_Drift
# from LIB_SIC_MODISAMSR2 import grab_SIC_MODISAMSR
# from LIB_access_ERA5 import shift_bylons