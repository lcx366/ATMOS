import numpy as np
from datetime import datetime,timedelta
import pandas as pd
from os import path,makedirs,remove
from pathlib import Path

from ..utils.try_download import wget_download

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
    
    swfile = direc + 'SW-All.csv'
    url = 'https://www.celestrak.com/SpaceData/SW-All.csv'

    if not path.exists(direc): makedirs(direc)
    if not path.exists(swfile):
        desc = "Downloading the Space Weather file '{:s}' from CELESTRAK".format('SW-All.csv')
        wget_download(url,swfile,desc)
    else:
        modified_time = datetime.fromtimestamp(path.getmtime(swfile))
        if datetime.now() > modified_time + timedelta(days=7):
            remove(swfile)
            desc = "Updating the Space Weather file '{:s}' from CELESTRAK".format('SW-All.csv')
            wget_download(url,swfile,desc)  
        else:
            print("The Space Weather file '{:s}' in {:s} is already the latest.".format('SW-All.csv',direc))   
    return swfile

def read_sw_nrlmsise00(swfile):
    '''
    Parse and read the space weather data

    Usage: 
    sw_obs_pre = read_sw_nrlmsise00(swfile)

    Inputs: 
    swfile -> [str] Path of the space weather data
    
    Outputs: 
    sw_obs_pre -> [2d str array] Content of the space weather data

    Examples:
    >>> swfile = 'sw-data/SW-All.csv'
    >>> sw_obs_pre = read_sw(swfile)
    >>> print(sw_obs_pre)
    [['2020' '01' '07' ... '72.4' '68.0' '71.0']
    ['2020' '01' '06' ... '72.4' '68.1' '70.9']
    ...
    ...
    ['1957' '10' '02' ... '253.3' '267.4' '231.7']
    ['1957' '10' '01' ... '269.3' '266.6' '230.9']]
    '''
    sw_df = pd.read_csv(swfile)  
    sw_df.dropna(subset=['C9'],inplace=True)
    # Sort from newest date to past
    sw_df.sort_values(by=['DATE'],ascending=False,inplace=True)
    sw_df.reset_index(drop=True,inplace=True)
    return sw_df

def get_sw(sw_df,t_ymd,hour):
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

    ymds = sw_df['DATE']
    j_, = np.where(sw_df['DATE'] == t_ymd)
    j = j_[0]
    f107A,f107,ap = sw_df.iloc[j]['F10.7_OBS_CENTER81'],sw_df.iloc[j+1]['F10.7_OBS'],sw_df.iloc[j]['AP_AVG']
    aph_tmp_b0 = sw_df.iloc[j]['AP1':'AP8']   
    i = int(np.floor_divide(hour,3))
    ap_c = aph_tmp_b0[i]
    aph_tmp_b1 = sw_df.iloc[j+1]['AP1':'AP8']
    aph_tmp_b2 = sw_df.iloc[j+2]['AP1':'AP8']
    aph_tmp_b3 = sw_df.iloc[j+3]['AP1':'AP8']
    aph_tmp = np.hstack((aph_tmp_b3,aph_tmp_b2,aph_tmp_b1,aph_tmp_b0))[::-1]
    apc_index = 7-i
    aph_c369 = aph_tmp[apc_index:apc_index+4]
    aph_1233 = np.average(aph_tmp[apc_index+4:apc_index+12])
    aph_3657 = np.average(aph_tmp[apc_index+12:apc_index+20])
    aph = np.hstack((ap,aph_c369,aph_1233,aph_3657))
    return f107A,f107,ap,aph