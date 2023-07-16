from astropy.utils import iers as iers_astropy

from .data_download import download_iers

def iers_load():

    # load the EOP file
    dir_iers,eop_file,leapsecond_file = download_iers() 
    iers_astropy.conf.auto_download = False
    iers_a = iers_astropy.IERS_A.open(eop_file)
    leapsecond = iers_astropy.LeapSeconds.from_iers_leap_seconds(leapsecond_file)
    eop_table = iers_astropy.earth_orientation_table.set(iers_a)