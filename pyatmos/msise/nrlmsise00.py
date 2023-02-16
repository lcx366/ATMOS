import numpy as np
from astropy.time import Time

from .nrlmsise00_subfunc import gtd7,gtd7d
from .spaceweather import get_sw
from ..utils.utils import wraplon,hms_conver
from ..class_atmos import ATMOS

def nrlmsise00(t,location,SW_OBS_PRE,aphmode=True):
    """
    NRLMSISE-00 is a semi-empirical, global reference atmospheric model of the Earth from ground level to the exosphere up to 2000 km.
    A primary use of this model is to aid predictions of satellite orbital decay due to the atmospheric drag. 

    Usage:
    nrl00 = nrlmsise00(t,(lat,lon,alt),swdata,aphmode = 'Aph')    

    Inputs:
    t -> [str] time(UTC)
    location -> [tuple/list] geodetic latitude[degree], longitude[degree], and altitude[km]

    Parameters:
    aphmode -> [bool, optional, default = True] whether to use the 3h geomagnetic index

    Output:
    nrl00 -> instance of class ATMOS, where its attributes include
        rho -> [float] total mass density[kg/m^3]
        T -> [tuple] local temperature[K]
        nd -> [dict] number density of components[1/m^3], including Helium(He), Oxygen(O), Oxygen(O2), Nitrogen(N),
        Nitrogen(N2), Argon(Ar), Hydrogen(H), Anomalous Oxygen(ANM O)

    Examples:
    >>> from pyatmos import download_sw_nrlmsise00,read_sw_nrlmsise00
    >>> # Download or update the space weather file from www.celestrak.com
    >>> swfile = download_sw_nrlmsise00() 
    >>> # Read the space weather data
    >>> swdata = read_sw_nrlmsise00(swfile) 
    >>>
    >>> from pyatmos import nrlmsise00
    >>> # Set a specific time and location
    >>> t = '2014-07-22 22:18:45' # time(UTC)   
    >>> lat,lon,alt = 25,102,600 # latitude, longitude in [degree], and altitude in [km]
    >>> nrl00 = nrlmsise00(t,(lat,lon,alt),swdata)
    >>> print(nrl00.rho) # [kg/m^3]
    >>> print(nrl00.T) # [K]
    >>> print(nrl00.nd) # composition in [1/m^3]     
    """

    lat,lon,h = location

    # calculate the altitude above sea level from height
    alt = h
        
    t = Time(t)
    t_ymd = t.isot.split('T')[0]
    t_yday = t.yday.split(':')
    year,doy = int(t_yday[0]),int(t_yday[1])
    hour,sec = hms_conver(int(t_yday[2]),int(t_yday[3]),float(t_yday[4]))
    lst = hour + wraplon(lon)/15
    if alt > 80:
        f107A,f107,ap,aph = get_sw(SW_OBS_PRE,t_ymd,hour)
    else:
        f107A,f107,ap,aph = 150,150,4,np.full(7,4)

    lon_wrap = wraplon(lon)
    inputp = {'doy':doy,'year':year,'sec':sec,'alt':alt,'g_lat':lat,'g_lon':lon_wrap,'lst':lst,\
              'f107A':f107A,'f107':f107,'ap':ap,'ap_a':aph}
    
    switches = np.ones(23)
    if aphmode: switches[8] = -1 # -1 indicates the use of 3h geomagnetic index
        
    if alt > 500:
        output = gtd7d(inputp,switches)
    else:
        output = gtd7(inputp,switches)

    inputp['g_lon'] = lon   
    params = {'Year':inputp['year'],'DOY':inputp['doy'],'SOD':inputp['sec'],'Lat':inputp['g_lat'],'Lon':inputp['g_lon'],'Alt':inputp['alt'],'LST':inputp['lst'],\
              'f107A':inputp['f107A'],'f107D':inputp['f107'],'ApD':inputp['ap'],'Ap3H':inputp['ap_a']}
    rho = output['d']['RHO']          
    T = (output['t']['TINF'],output['t']['TG'])
    nd = {'He':output['d']['He'],'O':output['d']['O'],'N2':output['d']['N2'],'O2':output['d']['O2'],'Ar':output['d']['AR'],'H':output['d']['H'],'N':output['d']['N'],'ANM O':output['d']['ANM O']}

    info = {'rho':rho,'T':output['t']['TG'],'nd':nd}

    return ATMOS(info)
