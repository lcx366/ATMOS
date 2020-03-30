# Welcome to ATMOS

The pyatmos package is an archive of scientific routines that aims to implement the estimation of atmospheric properties for various atmosphere models. Currently, feasible atmosphere models include:

1. International Standard Atmosphere(ISA) Model up to 86km
2. NRLMSISE-00

## How to install

pyatmos can be installed with

```sh
pip install pyatmos
```

## How to use

### International Standard Atmosphere

Calculate the ISA at an altitude(default is geometric) of 10km.

```python
>>> from pyatmos import isa
>>> isa(10)
{'temperature[K]': 223.25186489868483,
 'pressure[Pa]': 26499.756053713343,
 'density[kg/m^3]': 0.41350863360218376}
```

Calculate the ISA at a geopotential altitude of 50km.

```python
>>> isa(50,'geopotential')
{'temperature[K]': 270.65,
 'pressure[Pa]': 75.94476758456234,
 'density[kg/m^3]': 0.0009775244455727493}
```

Calculate the ISA at 90km.

```python
>>> isa(90)
Exception: geometric altitude should be in [-0.611, 86.0] km
>>> isa(90,'geopotential')    
Exception: geopotential altitude should be in [-0.610, 84.852] km    
```

### NRLMSISE-00

Get the space weather data

```python
>>> from pyatmos import download_sw,read_sw
>>> # Download or update the space weather file from www.celestrak.com
>>> swfile = download_sw() 
>>> # Read the space weather data
>>> swdata = read_sw(swfile) 
Updating the space weather data ... Finished
```

Calculate the temperatures, densities **not** including anomalous oxygen using the NRLMSISE-00 model at 70km, 25 degrees latitude, 102 degrees longitude on the date October 5, 2015 at 03:00:00 UTC.

```
>>> from pyatmos import nrlmsise00
>>> # Set a specific time and location
>>> t = '2015-10-05 03:00:00' # time(UTC)
>>> lat,lon = 25,102 # latitude and longitude [degree]
>>> alt = 70 # altitude [km]
>>> para_input,para_output = nrlmsise00(t,lat,lon,alt,swdata)
>>> print(para_input,'\n')
>>> print(para_output)
{'Year': 2015, 'DayOfYear': 278, 'SecondOfDay': 10800.0, 'Latitude[deg]': 25, 'Longitude[deg]': 102, 'Altitude[km]': 70, 'LocalSolarTime[hours]': 9.8, 'f107Average[10^-22 W/m^2/Hz]': 150, 'f107Daily[10^-22 W/m^2/Hz]': 150, 'ApDaily': 4, 'Ap3Hourly': array([4, 4, 4, 4, 4, 4, 4])} 

{'Density': {'He[1/m^3]': 9100292488300570.0, 'O[1/m^3]': 0, 'N2[1/m^3]': 1.3439413974205876e+21, 'O2[1/m^3]': 3.52551376755781e+20, 'AR[1/m^3]': 1.6044163757370681e+19, 'H[1/m^3]': 0, 'N[1/m^3]': 0, 'ANM O[1/m^3]': 0, 'RHO[kg/m^3]': 8.225931818480755e-05}, 'Temperature': {'TINF[K]': 1027.3184649, 'TG[K]': 219.9649472491653}}
```

Calculate the temperatures, densities **not** including anomalous oxygen using the NRLMSISE-00 model at 100km, -65 degrees latitude, -120 degrees longitude on the date July 8, 2004 at 10:30:50 UTC.

```
>>> t = '2004-07-08 10:30:50' 
>>> lat,lon,alt = -65,-120,100 
>>> para_input,para_output = nrlmsise00(t,lat,lon,alt,swdata)
>>> print(para_input,'\n')
>>> print(para_output)
{'Year': 2004, 'DayOfYear': 190, 'SecondOfDay': 37850.0, 'Latitude[deg]': -65, 'Longitude[deg]': -120, 'Altitude[km]': 100, 'LocalSolarTime[hours]': 2.5138888888888893, 'f107Average[10^-22 W/m^2/Hz]': 109.0, 'f107Daily[10^-22 W/m^2/Hz]': 79.3, 'ApDaily': 2, 'Ap3Hourly': array([2.   , 2.   , 2.   , 2.   , 2.   , 3.125, 4.625])} 

{'Density': {'He[1/m^3]': 119477307274636.89, 'O[1/m^3]': 4.1658304136233e+17, 'N2[1/m^3]': 7.521248904485598e+18, 'O2[1/m^3]': 1.7444969074975662e+18, 'AR[1/m^3]': 7.739495767665198e+16, 'H[1/m^3]': 22215754381448.5, 'N[1/m^3]': 152814261016.3964, 'ANM O[1/m^3]': 1.8278224834873257e-37, 'RHO[kg/m^3]': 4.584596293339505e-07}, 'Temperature': {'TINF[K]': 1027.3184649, 'TG[K]': 192.5868649143824}}
```

Calculate the temperatures, densities including anomalous oxygen using the NRLMSISE-00 model at 500km, 85 degrees latitude, 210 degrees longitude on the date February 15, 2010 at 12:18:37 UTC.

```
>>> t = '2010-02-15 12:18:37' 
>>> lat,lon,alt = 85,210,500 
>>> para_input,para_output = nrlmsise00(t,lat,lon,alt,swdata,omode='Oxygen')
>>> print(para_input,'\n')
>>> print(para_output)
{'Year': 2010, 'DayOfYear': 46, 'SecondOfDay': 44317.0, 'Latitude[deg]': 85, 'Longitude[deg]': 210, 'Altitude[km]': 500, 'LocalSolarTime[hours]': 2.310277777777779, 'f107Average[10^-22 W/m^2/Hz]': 83.4, 'f107Daily[10^-22 W/m^2/Hz]': 89.4, 'ApDaily': 14, 'Ap3Hourly': array([14.   ,  5.   ,  7.   ,  6.   , 15.   ,  5.375,  4.   ])} 

{'Density': {'He[1/m^3]': 2830075020953.2334, 'O[1/m^3]': 5866534735436.941, 'N2[1/m^3]': 59516979995.87239, 'O2[1/m^3]': 1558775273.2950978, 'AR[1/m^3]': 825564.7467165776, 'H[1/m^3]': 142697077779.00586, 'N[1/m^3]': 53473812381.891624, 'ANM O[1/m^3]': 4258921381.0652237, 'RHO[kg/m^3]': 1.790487924033088e-13}, 'Temperature': {'TINF[K]': 850.5598890315023, 'TG[K]': 850.5507885501303}}
```

Calculate the temperatures, densities including anomalous oxygen using the NRLMSISE-00 model at 900km, 3 degrees latitude, 5 degrees longitude on the date August 20, 2019 at 23:10:59 UTC. It uses not only Daily AP but also 3-hour AP magnetic index.

```
>>> t = '2019-08-20 23:10:59' 
>>> lat,lon,alt = 3,5,900 
>>> para_input,para_output = nrlmsise00(t,lat,lon,alt,swdata,omode='Oxygen',aphmode = 'Aph')
>>> print(para_input,'\n')
>>> print(para_output)
{'Year': 2019, 'DayOfYear': 232, 'SecondOfDay': 83459.0, 'Latitude[deg]': 3, 'Longitude[deg]': 5, 'Altitude[km]': 900, 'LocalSolarTime[hours]': 23.51638888888889, 'f107Average[10^-22 W/m^2/Hz]': 67.4, 'f107Daily[10^-22 W/m^2/Hz]': 67.7, 'ApDaily': 4, 'Ap3Hourly': array([4.   , 4.   , 3.   , 3.   , 5.   , 3.625, 3.5  ])} 

{'Density': {'He[1/m^3]': 74934329990.0412, 'O[1/m^3]': 71368139.39199762, 'N2[1/m^3]': 104.72048033793158, 'O2[1/m^3]': 0.09392848471935447, 'AR[1/m^3]': 1.3231114543012155e-07, 'H[1/m^3]': 207405192640.34592, 'N[1/m^3]': 3785341.821909535, 'ANM O[1/m^3]': 1794317839.638502, 'RHO[kg/m^3]': 8.914971667362366e-16}, 'Temperature': {'TINF[K]': 646.8157488121493, 'TG[K]': 646.8157488108872}}
```

## Change log
- **1.1.0 — Mar 29,  2020**
  - Added the International Standard Atmosphere(ISA) Model up to 86km 

## Next release

- Complete the help documentation
- Improve the code structure to make it easier to read
- Add other atmospheric models, such as the **U.S. Standard Atmosphere 1976(USSA1976)** or **Committee on Extension to the Standard Atmosphere(COESA)** up to 1000km, **Unofficial Australian Standard Atmosphere 2000(UASA2000)**, and the **Jacchia-Bowman 2008 Empirical Thermospheric Density Model(JB2008)**

## Reference

- U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, D.C. 
- [Public Domain Aeronautical Software](http://www.pdas.com/atmos.html) 
- https://gist.github.com/buzzerrookie/5b6438c603eabf13d07e
- https://ww2.mathworks.cn/help/aerotbx/ug/atmosisa.html

* [Original Fortran and C code](https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/)
* [MSISE-00 in Python and Matlab](https://github.com/space-physics/msise00)
* [NRLMSISE-00 Atmosphere Model - Matlab](https://ww2.mathworks.cn/matlabcentral/fileexchange/56253-nrlmsise-00-atmosphere-model?requestedDomain=zh)
* [NRLMSISE-00 Atmosphere Model - Aerospace Blockset](https://www.mathworks.com/help/aeroblks/nrlmsise00atmospheremodel.html?requestedDomain=)
* [NRLMSISE-00 Atmosphere Model - CCMC](https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php)
* [NRLMSISE-00 empirical model of the atmosphere: Statistical comparisons and scientific issues](http://onlinelibrary.wiley.com/doi/10.1029/2002JA009430/pdf)

