'''
pyatmos msise subpackage

This subpackage defines the following functions:

# nrlmsise00.py  - Implements the NRLMSISE 00 model

# spaceweather.py 

download_sw - Download or update the space weather file from www.celestrak.com

read_sw - Read the space weather file

get_sw - Extract the space weather data  

'''              
from .spaceweather import download_sw,read_sw