'''
pyatmos package

This package is an archive of scientific routines that implements the 
estimation of atmospheric properties for various atmosphere models, such
as exponential, coesa76, and nrlmsise00. 

The package mainly estimates density, temperature, pressure and other parameters 
of air at a set of specific altitudes. For atmosphere below 86 kilometers, it also 
calculates the speed of sound, viscosity, and thermal conductivity.
'''    

from .standardatmos.expo import expo
from .standardatmos.coesa76 import coesa76
from .msise.spaceweather import download_sw,read_sw
from .msise.nrlmsise00 import nrlmsise00