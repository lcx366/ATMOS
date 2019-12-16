# Welcome to ATMOS

The pyatmos package is an archive of scientific routines that can be used to handle atmospheric models. Currently, only nrlmsise00 is feasible, which is valid from altitude z = 0..1000 km.

## How to Install

```sh
pip install pyatmos
```

## Examples

```python
from pyatmos.msise import download_sw,read_sw
from pyatmos.atmosclasses import Coordinate

# Download or update the space weather file from www.celestrak.com
swfile = download_sw() 
# Read the space weather data
sw_obs_pre = read_sw(swfile) 

# Test 1
# Set a specific time and location
t = '2015-10-05 03:00:00' # time(UTC)
lat,lon = 25,102 # latitude and longitude [degree]
alt = 70 # altitude [km]
# Initialize a coordinate instance by a space-time point
st = Coordinate(t,lat,lon,alt)

para_input,para_output = st.nrlmsise00(sw_obs_pre)
print(para_input)
{'doy': 278, 'year': 2015, 'sec': 10800.0, 'alt': 70, 'g_lat': 25, 'g_long': 102, 'lst': 9.8, 'f107A': 150, 'f107': 150, 'ap': 4, 'ap_a': array([4, 4, 4, 4, 4, 4, 4])}
print(para_output)
{'d': {'He': 9100292488300570.0, 'O': 0, 'N2': 1.3439413974205876e+21, 'O2': 3.52551376755781e+20, 'AR': 1.6044163757370681e+19, 'RHO': 8.225931818480755e-05, 'H': 0, 'N': 0, 'ANM O': 0}, 't': {'TINF': 1027.3184649, 'TG': 219.9649472491653}}

# Test 2
t = '2004-07-08 10:30:50' 
lat,lon,alt = -65,-120,100 
st = Coordinate(t,lat,lon,alt)
para_input,para_output = st.nrlmsise00(sw_obs_pre)
print(para_input)
{'doy': 190, 'year': 2004, 'sec': 37850.0, 'alt': 100, 'g_lat': -65, 'g_long': -120, 'lst': 2.5138888888888893, 'f107A': 109.0, 'f107': 79.3, 'ap': 2, 'ap_a': array([2.   , 2.   , 2.   , 2.   , 2.   , 3.125, 4.625])}
print(para_output)
{'d': {'He': 119477307274636.89, 'O': 4.1658304136233e+17, 'N2': 7.521248904485598e+18, 'O2': 1.7444969074975662e+18, 'AR': 7.739495767665198e+16, 'RHO': 4.584596293339505e-07, 'H': 22215754381448.5, 'N': 152814261016.3964, 'ANM O': 1.8278224834873257e-37}, 't': {'TINF': 1027.3184649, 'TG': 192.5868649143824}}

# Test 3
t = '2010-02-15 12:18:37' 
lat,lon,alt = 85,210,500 
st = Coordinate(t,lat,lon,alt)
para_input,para_output = st.nrlmsise00(sw_obs_pre,'NoOxygen','Aph')
print(para_input)
{'doy': 46, 'year': 2010, 'sec': 44317.0, 'alt': 500, 'g_lat': 85, 'g_long': 210, 'lst': 2.310277777777779, 'f107A': 83.4, 'f107': 89.4, 'ap': 14, 'ap_a': array([14.   ,  5.   ,  7.   ,  6.   , 15.   ,  5.375,  4.   ])}
print(para_output)
{'d': {'He': 3314507585382.5425, 'O': 3855595951659.0874, 'N2': 19285497858.028534, 'O2': 395599656.3119481, 'AR': 146073.85956102316, 'RHO': 1.2650700238089615e-13, 'H': 171775437382.8238, 'N': 38359828672.39737, 'ANM O': 5345258193.554493}, 't': {'TINF': 776.3155804924045, 'TG': 776.3139192714452}}

# Test 4
t = '2019-08-20 23:10:59' 
lat,lon,alt = 3,5,900 
st = Coordinate(t,lat,lon,alt)
para_input,para_output = st.nrlmsise00(sw_obs_pre,aphmode = 'Aph')
print(para_input)
{'doy': 232, 'year': 2019, 'sec': 83459.0, 'alt': 900, 'g_lat': 3, 'g_long': 5, 'lst': 23.51638888888889, 'f107A': 67.4, 'f107': 67.7, 'ap': 4, 'ap_a': array([4.   , 4.   , 3.   , 3.   , 5.   , 3.625, 3.5  ])}
print(para_output)
{'d': {'He': 74934329990.0412, 'O': 71368139.39199762, 'N2': 104.72048033793158, 'O2': 0.09392848471935447, 'AR': 1.3231114543012155e-07, 'RHO': 8.914971667362366e-16, 'H': 207405192640.34592, 'N': 3785341.821909535, 'ANM O': 1794317839.638502}, 't': {'TINF': 646.8157488121493, 'TG': 646.8157488108872}}
```

## Reference

* [Original Fortran and C code](https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/)
* [MSISE-00 in Python and Matlab](https://github.com/space-physics/msise00)
* [NRLMSISE-00 Atmosphere Model - Matlab](https://ww2.mathworks.cn/matlabcentral/fileexchange/56253-nrlmsise-00-atmosphere-model?requestedDomain=zh)
* [NRLMSISE-00 Atmosphere Model - Aerospace Blockset](https://www.mathworks.com/help/aeroblks/nrlmsise00atmospheremodel.html?requestedDomain=)
* [NRLMSISE-00 Atmosphere Model - CCMC](https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php)
* [NRLMSISE-00 empirical model of the atmosphere: Statistical comparisons and scientific issues](http://onlinelibrary.wiley.com/doi/10.1029/2002JA009430/pdf)

