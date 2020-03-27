'''
pyatmos msise subpackage

This subpackage defines the following functions:

# ==================== nrlmisise-00  model =================== #

nrlmsise00  - Implements the NRLMSISE 00 model

# =================== spaceweather functions ================= #

download_sw - Download or update the space weather file from www.celestrak.com

read_sw     - Read the space weather file

get_sw      - Extract the space weather data  

# ===================== utility functions ==================== #

wraplon     - Wrap a longitude in range of [0,360] to [-180,180]

hms2s       - Convert hour/minute/second to seconds

hms2h       - Convert hour/minute/second to hours  

'''              
from .spaceweather import download_sw,read_sw