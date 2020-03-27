'''
pyatmos package

This package is an archive of scientific routines that implements the estimation of atmospheric properties for various atmosphere models.  
'''    

from .StandardAtmosphere.StandardAtmosphere import isa
from .msise.spaceweather import download_sw,read_sw
from .msise.nrlmsise00 import nrlmsise00