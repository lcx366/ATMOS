import numpy as np
from datetime import datetime,timedelta
from os import path,makedirs,remove
from pathlib import Path

from ..utils.try_download import tqdm_request

def download_sw_nrlmsise00(direc=None):
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
    
    swfile = direc + 'SW-All.txt'
    url = 'https://www.celestrak.com/SpaceData/SW-All.txt'

    if not path.exists(direc): makedirs(direc)
    if not path.exists(swfile):
        desc = 'Downloading the space weather data {:s} from CELESTRAK'.format('SW-All.txt')
        tqdm_request(url,direc,'SW-All.txt',desc)
    else:
        modified_time = datetime.fromtimestamp(path.getmtime(swfile))
        if datetime.now() > modified_time + timedelta(days=1):
            remove(swfile)
            desc = 'Updating the space weather data {:s} from CELESTRAK'.format('SW-All.txt')
            tqdm_request(url,direc,'SW-All.txt',desc)   
        else:
            print('The space weather data in {:s} is already the latest.'.format(direc))   
    return swfile
 
def read_sw_nrlmsise00(swfile):
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
    sw_data = open(swfile,'r').readlines()
    SW_OBS,SW_PRE = [],[]
    flag1 = flag2 = 0
    for line in sw_data:
        if line.startswith('BEGIN OBSERVED'): 
            flag1 = 1
            continue
        if line.startswith('END OBSERVED'): flag1 = 0 
        if flag1 == 1: 
            sw_p = line.split()
            if len(sw_p) == 30:
                del sw_p[24]
            elif len(sw_p) == 31: 
                sw_p = np.delete(sw_p,[23,25]) 
            else: 
                sw_p = np.delete(sw_p,[23,24,25,27])
            SW_OBS.append(sw_p)
            
        if line.startswith('BEGIN DAILY_PREDICTED'): 
            flag2 = 1
            continue    
        if line.startswith('END DAILY_PREDICTED'): break 
        if flag2 == 1: SW_PRE.append(line.split())    
    SW_OBS_PRE = np.vstack((np.array(SW_OBS),np.array(SW_PRE)))   
    # inverse sort
    SW_OBS_PRE = np.flip(SW_OBS_PRE,0).astype(dtype='<U8')
    ymds = np.apply_along_axis(''.join, 1, SW_OBS_PRE[:,:3])
    SW_OBS_PRE = np.insert(SW_OBS_PRE[:,3:],0,ymds,axis=1)
    return SW_OBS_PRE 
 
def get_sw(SW_OBS_PRE,t_ymd,hour):
    '''
    Extract the necessary parameters describing the solar activity and geomagnetic activity from the space weather data.

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

    ymds = SW_OBS_PRE[:,0]
    j_, = np.where(''.join(t_ymd) == ymds)
    j = j_[0]
    f107A,f107,ap = float(SW_OBS_PRE[j,25]),float(SW_OBS_PRE[j+1,24]),int(SW_OBS_PRE[j,20])
    aph_tmp_b0 = SW_OBS_PRE[j,12:20]   
    i = int(np.floor_divide(hour,3))
    ap_c = aph_tmp_b0[i]
    aph_tmp_b1 = SW_OBS_PRE[j+1,12:20]
    aph_tmp_b2 = SW_OBS_PRE[j+2,12:20]
    aph_tmp_b3 = SW_OBS_PRE[j+3,12:20]
    aph_tmp = np.hstack((aph_tmp_b3,aph_tmp_b2,aph_tmp_b1,aph_tmp_b0))[::-1].astype(np.float)
    apc_index = 7-i
    aph_c369 = aph_tmp[apc_index:apc_index+4]
    aph_1233 = np.average(aph_tmp[apc_index+4:apc_index+12])
    aph_3657 = np.average(aph_tmp[apc_index+12:apc_index+20])
    aph = np.hstack((ap,aph_c369,aph_1233,aph_3657))
    return f107A,f107,ap,aph
