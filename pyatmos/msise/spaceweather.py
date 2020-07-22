# -------------------------------------------------------------------- #
# ------------------------  space weather  --------------------------- #
# -------------------------------------------------------------------- #

'''
pyatmos spaceweather

This submodule defines the following functions:

download_sw - Download or update the space weather file from www.celestrak.com

read_sw     - Read the space weather file

get_sw      - Extract the space weather data    
'''

import numpy as np
from datetime import datetime,timedelta
from os import path,makedirs,remove
from pathlib import Path
from urllib.request import urlretrieve
from colorama import Fore
from .tqdmupto import TqdmUpTo

# =================== download and update sw data  =================== #

def download_sw(direc=None):
    '''
    Download or update the space weather file from www.celestrak.com

    Usage: 
    swfile = download_sw([direc])

    Inputs: 
    direc -> [str, optionanl, default = $HOME+'/src/sw-data/'] Directory for storing sw file
    
    Outputs: 
    swfile -> [str] Path of sw file

    Examples:
    >>> swfile = download_sw()
    Downloading the latest space weather data ... Finished
    >>> print(swfile)
    /Users/lichunxiao/src/sw-data/SW-All.txt
    >>> swfile = download_sw('sw-data/')
    Downloading the latest space weather data ... Finished
    >>> swfile = download_sw('sw-data/')
    The existing space weather data is already up to date
    >>> print(swfile)
    sw-data/SW-All.txt
    '''
    
    if direc is None:
        home = str(Path.home())
        direc = home + '/src/sw-data/'
    
    swfile = direc + 'SW-All.txt'
    url = 'https://www.celestrak.com/SpaceData/SW-All.txt'

    if not path.exists(direc): makedirs(direc)
    bar_format = "{l_bar}%s{bar}%s{r_bar}" % (Fore.BLUE, Fore.RESET)
    if not path.exists(swfile):
        desc = 'Downloading the latest space weather data from CELESTRAK'
        with TqdmUpTo(unit='B', unit_scale=True, desc=desc,bar_format = bar_format) as t:
            urlretrieve(url, swfile,reporthook=t.update_to)
            t.total = t.n

    else:
        modified_time = datetime.fromtimestamp(path.getmtime(swfile))
        if datetime.now() > modified_time + timedelta(days=1):
            remove(swfile)
            desc = 'Updating the space weather data from CELESTRAK'
            with TqdmUpTo(unit='B', unit_scale=True, desc=desc,bar_format = bar_format) as t:
                urlretrieve(url, swfile,reporthook=t.update_to)    
                t.total = t.n
        else:
            print('The space weather data in {:s} is already the latest.'.format(direc))   
    return swfile

# =========================== read sw file ========================== #

def read_sw(swfile):
    '''
    Parse and read the space weather file

    Usage: 
    sw_obs_pre = read_sw(swfile)

    Inputs: 
    swfile -> [str] Path of sw file
    
    Outputs: 
    sw_obs_pre -> [2d str array] sw data

    Examples:
    >>> swfile = 'sw-data/SW-All.txt'
    >>> sw_obs_pre = read_sw(swfile)
    >>> print(sw_obs_pre)
    [['2020' '01' '07' ... '72.4' '68.0' '71.0']
    ['2020' '01' '06' ... '72.4' '68.1' '70.9']
    ['2020' '01' '05' ... '72.4' '68.2' '70.9']
    ...
    ['1957' '10' '03' ... '266.3' '268.1' '232.7']
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
    SW_OBS_PRE = np.flip(SW_OBS_PRE,0)
    return SW_OBS_PRE

# ========================== extract sw data ========================== #   

def get_sw(SW_OBS_PRE,t_ymd,hour):
    j = 0
    for ymd in SW_OBS_PRE[:,:3]:
        if np.array_equal(t_ymd,ymd): break
        j+=1 
    f107A,f107,ap = float(SW_OBS_PRE[j,27]),float(SW_OBS_PRE[j+1,26]),int(SW_OBS_PRE[j,22])
    aph_tmp_b0 = SW_OBS_PRE[j,14:22]   
    i = int(np.floor_divide(hour,3))
    ap_c = aph_tmp_b0[i]
    aph_tmp_b1 = SW_OBS_PRE[j+1,14:22]
    aph_tmp_b2 = SW_OBS_PRE[j+2,14:22]
    aph_tmp_b3 = SW_OBS_PRE[j+3,14:22]
    aph_tmp = np.hstack((aph_tmp_b3,aph_tmp_b2,aph_tmp_b1,aph_tmp_b0))[::-1].astype(np.float)
    apc_index = 7-i
    aph_c369 = aph_tmp[apc_index:apc_index+4]
    aph_1233 = np.average(aph_tmp[apc_index+4:apc_index+12])
    aph_3657 = np.average(aph_tmp[apc_index+12:apc_index+20])
    aph = np.hstack((ap,aph_c369,aph_1233,aph_3657))
    return f107A,f107,ap,aph
