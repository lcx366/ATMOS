import numpy as np
from datetime import datetime,timedelta
from os import path,makedirs,remove
from pathlib import Path

from ..utils.try_download import wget_download
from ..utils import Const

def download_sw_jb2008(direc=None):
    '''
    Download or update the space weather data from www.celestrak.com

    Usage: 
    swfile = download_sw([direc])

    Inputs: 
    direc -> [str, optional] Directory for storing the space weather data
    
    Outputs: 
    swfile -> [str] Path of the space weather data

    Examples:
    >>> swfile = download_sw()
    >>> swfile = download_sw('sw-data/')
    '''
    
    if direc is None:
        home = str(Path.home())
        direc = home + '/src/sw-data/'
    
    swfile1 = direc + 'SOLFSMY.TXT'
    url1 = 'https://sol.spacenvironment.net/JB2008/indices/SOLFSMY.TXT'

    # The file DTCFILE.TXT is generated from DSTFILE.TXT and SOLRESAP.TXT
    swfile2 = direc + 'DTCFILE.TXT'
    url2 = 'https://sol.spacenvironment.net/JB2008/indices/DTCFILE.TXT'

    if not path.exists(direc): makedirs(direc)
    if not path.exists(swfile1):
        desc1 = "Downloading the Space Weather file '{:s}' from Space Environment Technologies(SET)".format('SOLFSMY.TXT')
        desc2 = "Downloading the Space Weather file '{:s}' from Space Environment Technologies(SET)".format('DTCFILE.TXT')
        wget_download(url1,swfile1,desc1)
        wget_download(url2,swfile2,desc2)
    else:
        modified_time = datetime.fromtimestamp(path.getmtime(swfile1))
        if datetime.now() > modified_time + timedelta(days=1):
            remove(swfile1)
            remove(swfile2)
            desc1 = "Updating the Space weather data '{:s}' from Space Environment Technologies(SET)".format('SOLFSMY.TXT')
            desc2 = "Updating the Space weather data '{:s}' from Space Environment Technologies(SET)".format('DTCFILE.TXT')
            wget_download(url1,swfile1,desc1)
            wget_download(url2,swfile2,desc2)
        else:
            print("The Space Weather files '{:s}' and '{:s}' in {:s} are already the latest.".format('SOLFSMY.TXT','DTCFILE.TXT',direc))   
    return [swfile1,swfile2]
 
def read_sw_jb2008(swfile):
    '''
    Parse and read the space weather data

    Usage: 
    sw_obs_pre = read_sw(swfile)

    Inputs: 
    swfile -> [str] Path of the space weather data
    
    Outputs: 
    sw_obs_pre -> [2d str array] Content of the space weather data

    Examples:
    >>> swfile = 'sw-data/SW-All.txt'
    >>> sw_obs_pre = read_sw(swfile)
    >>> print(sw_obs_pre)
    [['2020' '01' '07' ... '72.4' '68.0' '71.0']
    ['2020' '01' '06' ... '72.4' '68.1' '70.9']
    ...
    ...
    ['1957' '10' '02' ... '253.3' '267.4' '231.7']
    ['1957' '10' '01' ... '269.3' '266.6' '230.9']]
    '''
    swfile1,swfile2 = swfile
    sw_data1 = np.loadtxt(swfile1,usecols=range(2,11))
    sw_data2 = np.loadtxt(swfile2,usecols=range(3,27),dtype=int)

    return (sw_data1,sw_data2)
 
def get_sw(sw_data,t_mjd):
    '''
    Extract the necessary parameters describing the solar activity and geomagnetic activity from the space weather data.
    INPUT T1950 AND READ FILE FOR VALUES FOR F10, S10, M10, AND Y10
    READ ONE TIME AND STORE ALL FILE VALUES IN COMMON FOR RETRIEVAL
    Usage: 
    f107A,f107,ap,aph = get_sw(SW_OBS_PRE,t_ymd,hour)

    Inputs: 
    SW_OBS_PRE -> [2d str array] Content of the space weather data
    t_ymd -> [str array or list] ['year','month','day']
    hour -> []
    
    Outputs: 
    f107A -> [float] 81-day average of F10.7 flux
    f107 -> [float] daily F10.7 flux for previous day
    ap -> [int] daily magnetic index 
    aph -> [float array] 3-hour magnetic index 

    Examples:
    >>> f107A,f107,ap,aph = get_sw(SW_OBS_PRE,t_ymd,hour)
    '''
    sw_data1,sw_data2 = sw_data
    sw_mjd = sw_data1[:,0] - 2400000.5
    J_, = np.where(sw_mjd-0.5 < t_mjd)
    j = J_[-1]

    # USE 1 DAY LAG FOR F10 AND S10 FOR JB2008
    dlag = j-1
    F10,F10B,S10,S10B = sw_data1[dlag,1:5] 

    # USE 2 DAY LAG FOR M10 FOR JB2008
    dlag = j-2
    M10,M10B = sw_data1[dlag,5:7] 

    # USE 5 DAY LAG FOR Y10 FOR JB2008
    dlag = j-5
    Y10,Y10B = sw_data1[dlag,7:9] 
    
    t_dmjd = t_mjd - sw_mjd[j] + 0.5
    x = Const.x
    y = sw_data2[j:j+2].flatten()
    DTCVAL = np.interp(t_dmjd,x,y)  
        
    return F10,F10B,S10,S10B,M10,M10B,Y10,Y10B,DTCVAL
