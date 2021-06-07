import numpy as np
from astropy.coordinates import get_sun
from astropy.time import Time

from .JB2008_subfunc import JB2008
from .spaceweather import get_sw
from ..utils.utils import ydhms_days
from ..class_atmos import ATMOS

def jb2008(t,location,swdata):
    """
    JB2008 is an empirical, global reference atmospheric model of the Earth from 90 km above the sea level to the exosphere up to 2500 km.
    A primary use of this model is to aid predictions of satellite orbital decay due to the atmospheric drag. 

    Usage:
    jb08 = jb2008(t,(lat,lon,alt),swdata)   

    Inputs:
    t -> [str] time(UTC)
    location -> [tuple/list] geodetic latitude[degree], longitude[degree], and altitude[km]

    Output:
    jb08 -> instance of class ATMOS, where its attributes include
        rho -> [float] total mass density[kg/m^3]
        T -> [tuple] local temperature[K]

    Examples:
    >>> from pyatmos import download_sw_jb2008,read_sw_jb2008
    >>> # Download or update the space weather file from https://sol.spacenvironment.net
    >>> swfile = download_sw_jb2008() 
    >>> # Read the space weather data
    >>> swdata = read_sw_jb2008(swfile)  
    >>>
    >>> from pyatmos import jb2008
    >>> # Set a specific time and location
    >>> t = '2014-07-22 22:18:45' # time(UTC) 
    >>> lat,lon,alt = 25,102,600 # latitude, longitude in [degree], and altitude in [km]
    >>> jb08 = jb2008(t,(lat,lon,alt),swdata)
    >>> print(jb08.rho) # [kg/m^3]
    >>> print(jb08.T) # [K]     
    """

    lat,lon,h = location
    t = Time(t,location=(str(lon)+'d',str(lat)+'d'))
    AMJD = t.mjd
    ydhms = np.array(t.yday.split(':'),dtype=float)
    YRDAY = ydhms_days(ydhms)
    sunpos = get_sun(t) 
    SUN = (sunpos.ra.rad,sunpos.dec.rad)
    sat_ra = t.sidereal_time('mean').rad
    SAT = (sat_ra,np.deg2rad(lat),h)

    F10,F10B,S10,S10B,M10,M10B,Y10,Y10B,DTCVAL = get_sw(swdata,AMJD)
    TEMP,RHO = JB2008(AMJD,YRDAY,SUN,SAT,F10,F10B,S10,S10B,M10,M10B,Y10,Y10B,DTCVAL)

    info = {'rho':RHO,'T':TEMP[1]}

    return ATMOS(info)
