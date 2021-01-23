# Welcome to ATMOS

[![PyPI version shields.io](https://img.shields.io/pypi/v/pyatmos.svg)](https://pypi.python.org/pypi/pyatmos/) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/pyatmos.svg)](https://pypi.python.org/pypi/pyatmos/) [![PyPI status](https://img.shields.io/pypi/status/pyatmos.svg)](https://pypi.python.org/pypi/pyatmos/) [![GitHub contributors](https://img.shields.io/github/contributors/lcx366/ATMOS.svg)](https://GitHub.com/lcx366/ATMOS/graphs/contributors/) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/lcx366/ATMOS/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/lcx366/ATMOS.svg)](https://github.com/lcx366/ATMOS/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/pystmos/badge/?version=latest)](http://pyatmos.readthedocs.io/?badge=latest) [![Build Status](https://travis-ci.org/lcx366/ATMOS.svg?branch=master)](https://travis-ci.org/lcx366/ATMOS)

This package is an archive of scientific routines that implements the estimation of atmospheric properties for various atmosphere models, such as Exponential, COESA76, and NRLMSISE-00. The package mainly estimates density, temperature, and pressure of air at a set of specific altitudes.

## How to install

On Linux, macOS and Windows architectures, the binary wheels can be installed using pip by executing one of the following commands:

```sh
pip install pyatmos
pip install pyatmos --upgrade # to upgrade a pre-existing installation
```

## How to use

#### Exponential

```python
>>> from pyatmos import expo
>>> rhos_geom = expo([0,20,40,60,80]) # geometric altitudes by default
>>> print(rhos_geom) # [kg/m^3]
>>> rhos_geop = expo([0,20,40,60,80],'geopotential') # geopotential altitudes
>>> print(rhos_geop)

[1.22500000e+00 7.76098911e-02 3.97200000e-03 3.20600000e-04
 1.90500000e-05]
[1.22500000e+00 7.69385063e-02 3.84131212e-03 2.97747719e-04
 1.59847603e-05]
```

#### COESA 1976

```python
>>> from pyatmos import coesa76
>>> rhos_geom,Ts_geom,Ps_geom = coesa76([0,20,40,60,80]) 
>>> print(rhos_geom) # [kg/m^3]
>>> rhos_geop,Ts_geop,Ps_geop = coesa76([0,20,40,60,80],'geopotential')
>>> print(rhos_geop) # [kg/m^3]
>>> rhos_geom,Ts_geom,Ps_geom = coesa76([100,300,500,700,900]) 
>>> print(rhos_geom) # [kg/m^3]

[1.22499916e+00 8.89079563e-02 3.99535051e-03 3.09628985e-04
 1.84514759e-05]
[1.22499916e+00 8.80348036e-02 3.85100688e-03 2.88320680e-04
 1.57005388e-05]
[5.60184300e-07 1.91512264e-11 5.21285933e-13 3.06944380e-14
 5.75807856e-15]  
```

#### NRLMSISE-00

*Before using NRLMSISE-00, the space weather data needs to be prepared in advance.*

```python
>>> from pyatmos import download_sw,read_sw
>>> # Download or update the space weather file from www.celestrak.com
>>> swfile = download_sw() 
>>> # Read the space weather data
>>> swdata = read_sw(swfile) 
```

Calculate the temperature, density at [25N, 102E, 20km] at 03:00:00 UTC on October 5, 2015 with anomalous oxygen and 3h-geomagnetic index.

```
>>> from pyatmos import nrlmsise00
>>> # Set a specific time and location
>>> t = '2015-10-05 03:00:00' # time(UTC) 
>>> lat,lon,alt = 25,102,600 # latitude, longitude in [degree], and altitude in [km]
>>> params,rho,T,nd = nrlmsise00(t,(lat,lon,alt),swdata) # aphmode=True
>>> print(params)
>>> print(rho) 
>>> print(T) 
>>> print(nd)

{'Year': 2015, 'DOY': 278, 'SOD': 10800.0, 'Lat': 25, 'Lon': 102, 'Alt': 600, 'LST': 9.8, 'f107A': 104.4, 'f107D': 82.6, 'ApD': 18, 'Ap3H': array([18.   , 22.   , 22.   , 22.   ,  7.   , 15.25 ,  9.375])}
6.416602651204796e-14
(853.466244160143, 853.4647165799171)
{'He': 2388916051039.6826, 'O': 1758109067905.8027, 'N2': 2866987110.5606275, 'O2': 22411077.605527952, 'Ar': 4351.013995142538, 'H': 155026672753.3203, 'N': 46719306249.863495, 'ANM O': 4920851253.780525}
```

**Note: The range of longitude is [0,360] by default, and the west longitude can also be expressed as a negative number.**

## Change log
- **1.2.1 — Jan 22, 2021**
  - Added **Exponential Atmosphere** up to 1000 km
  - Added **Committee on Extension to the Standard Atmosphere(COESA)** up to 1000 km
  - Completed part of the help documentation for NRLMSISE-00
  - Improved the code structure to make it easier to read
- **1.1.2 — Jul 26, 2020**
  - Added colored-progress bar for downloading data
- **1.1.0 — Mar 29,  2020**
  - Added the International Standard Atmosphere(ISA) Model up to 86kms  

## Next release

- Complete the help documentation for NRLMSISE-00
- Add other atmospheric models, such as the **Earth Global Reference Atmospheric Model(Earth-GRAM) 2016**, and the **Jacchia-Bowman 2008 Empirical Thermospheric Density Model(JB2008)**

## Reference

- U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, D.C. 
- [Public Domain Aeronautical Software](http://www.pdas.com/atmos.html) 
- https://gist.github.com/buzzerrookie/5b6438c603eabf13d07e
- https://ww2.mathworks.cn/help/aerotbx/ug/atmosisa.html

- [Original Fortran and C code](https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/)
- [MSISE-00 in Python and Matlab](https://github.com/space-physics/msise00)
- [NRLMSISE-00 Atmosphere Model - Matlab](https://ww2.mathworks.cn/matlabcentral/fileexchange/56253-nrlmsise-00-atmosphere-model?requestedDomain=zh)
- [NRLMSISE-00 Atmosphere Model - Aerospace Blockset](https://www.mathworks.com/help/aeroblks/nrlmsise00atmospheremodel.html?requestedDomain=)
- [NRLMSISE-00 Atmosphere Model - CCMC](https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php)
- [NRLMSISE-00 empirical model of the atmosphere: Statistical comparisons and scientific issues](http://onlinelibrary.wiley.com/doi/10.1029/2002JA009430/pdf)
- [ATMOSPHERIC MODELS](http://www.braeunig.us/space/atmmodel.htm)
- [poliastro-Atmosphere module](https://docs.poliastro.space/en/stable/api/safe/atmosphere/atmosphere_index.html?highlight=nrlmsise#famous-atmospheric-models)

