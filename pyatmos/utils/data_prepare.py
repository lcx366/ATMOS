from astropy.utils import iers as iers_astropy

from .data_download import download_iers

def iers_load(enable_iers_load=True):
    """
    Loads the EOP (Earth Orientation Parameters) and Leap Second files from IERS.
    This function downloads the necessary files if they are not found locally, and then sets up the Astropy libraries to use this data.
    """
    if enable_iers_load:
        # Download the EOP and Leap Second files
        dir_iers,eop_file,leapsecond_file = download_iers()
        # Load IERS data for Astropy
        iers_astropy.conf.auto_download = False # Prevent automatic IERS download by Astropy
        iers_a = iers_astropy.IERS_A.open(eop_file) # Load IERS data
        leapsecond = iers_astropy.LeapSeconds.from_iers_leap_seconds(leapsecond_file) # Load Leap Second data
        eop_table = iers_astropy.earth_orientation_table.set(iers_a) # Configure Astropy to use IERS data
