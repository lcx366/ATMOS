from datetime import datetime,timedelta
from os import path,makedirs,remove
from pathlib import Path

from .try_download import wget_download 

def download_iers(out_days=7,dir_to=None):
    """
    Download or update the Earth Orientation Parameters(EOP) file and Leap Second file from IERS

    Usage: 
        >>> dir_to,dir_eop_file,dir_leapsecond_file = download_iers()
    Inputs: 
        out_days -> [int, optional, default = 7] Updating cycle of the IERS files
        dir_to   -> [str, optional, default = None] Directory for storing EOP file
    Outputs: 
        dir_to -> [str] Directory of the IERS files
        dir_eop_file -> [str] Path of the EOP file
        dir_leapsecond_file -> [str] Path of the Leap Second file
    """
    if dir_to is None:
        home = str(Path.home())
        dir_to = home + '/src/iers/'
    
    eop_file = 'finals2000A.all'
    leapsecond_file = 'Leap_Second.dat'
    dir_eop_file = dir_to + eop_file
    dir_leapsecond_file = dir_to + leapsecond_file

    url_eop = 'https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all'
    url_leapsecond = 'https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat'

    if not path.exists(dir_to): makedirs(dir_to)
    if not path.exists(dir_eop_file):
        desc = "Downloading the latest EOP file '{:s}' from IERS".format(eop_file)
        wget_out = wget_download(url_eop,dir_eop_file,desc)
    else:
        modified_time = datetime.fromtimestamp(path.getmtime(dir_eop_file))
        if datetime.now() > modified_time + timedelta(days=out_days):
            remove(dir_eop_file)
            desc = "Updating the EOP file '{:s}' from IERS".format(eop_file)
            wget_out = wget_download(url_eop,dir_eop_file,desc)
        else:
            print("The EOP file '{:s}' in {:s} is already the latest.".format(eop_file,dir_to)) 

    if not path.exists(dir_leapsecond_file):
        desc = "Downloading the latest Leap Second file '{:s}' from IERS".format(leapsecond_file)
        wget_out = wget_download(url_leapsecond,dir_leapsecond_file,desc)
    else:
        modified_time = datetime.fromtimestamp(path.getmtime(dir_leapsecond_file))
        if datetime.now() > modified_time + timedelta(days=out_days):
            remove(dir_leapsecond_file)
            desc = "Updating the Leap Second file '{:s}' from IERS".format(leapsecond_file)
            wget_out = wget_download(url_leapsecond,dir_leapsecond_file,desc)
        else:
            print("The Leap Second file '{:s}' in {:s} is already the latest.".format(leapsecond_file,dir_to))        

    return dir_to,dir_eop_file,dir_leapsecond_file    