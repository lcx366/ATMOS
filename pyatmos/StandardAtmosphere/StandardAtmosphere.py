import numpy as np

# physical constants:
g = 9.80665 # standard gravity acceleration in [m/s^2]
R = 8314.32/28.9644 # gas constant for dry air in [J/kg/K] from the 1976 United States Standard Atmosphere(USSA) Model
R_e = 6371000.8 # Earth's volumetric radius in [m]
# temperature and pressure at the Mean Sea Level(MSL)
t_msl = 288.15 # temperature in [K]
p_msl = 101325 # pressure in [Pa]
h_msl = 0 # altitude in [m]
    
def lapse_tp(t_lower, p_lower, lr, h_lower, h_upper):
    '''
    Calculate the temperature and pressure at a given geopotential altitude above base of a specific layer.
    The temperature is computed by the linear interpolation with the slope defined by lapse rates. 
    The pressure is computed from the hydrostatic equations and the perfect gas law. 
    See the detailed documentation in http://www.pdas.com/hydro.pdf

    Usage:
    [t_upper, p_upper] = lapse_tp(t_lower, p_lower, lr, h_lower, h_upper)

    Inputs:
    t_lower -> [float] temperature in [K] at the lower boundary of the subset in the specific layer
    p_lower -> [float] pressure in [Pa] at the lower boundary of the subset in the specific layer
    lr -> [float] lapse rate in [K/m] for the specific layer
    h_lower -> [float] geopotential altitude in [m] at the lower boundary of the subset in the specific layer
    h_upper -> [float] geopotential altitude in [m] at the upper boundary of the subset in the specific layer
    
    Outputs:
    t1 -> [float] temperature in [K] at the upper boundary of the subset in the specific layer
    p1 -> [float] pressure in [Pa] at the upper boundary of the subset in the specific layer
    
    Reference: Public Domain Aeronautical Software(http://www.pdas.com/atmos.html) 
               https://gist.github.com/buzzerrookie/5b6438c603eabf13d07e
    '''
    if lr == 0:
        t_upper = t_lower
        p_upper = p_lower * np.exp(-g / R / t_lower * (h_upper - h_lower))
    else:
        t_upper = t_lower + lr * (h_upper - h_lower)
        p_upper = p_lower * (t_upper / t_lower) ** (-g / lr / R)

    return t_upper,p_upper

def isa(altitude,altitude_type='geometric'):
    '''
    Implements the International Standard Atmosphere(ISA) Model up to 86km. 
    The standard atmosphere is defined as a set of layers by specified geopotential altitudes and lapse rates.
    The temperature is computed by linear interpolation with the slope defined by the lapse rate. 
    The pressure is computed from the hydrostatic equations and the perfect gas law; the density follows from the perfect gas law. 

    Usage:
    [temperature, pressure, density] = isa(altitude)
    [temperature, pressure, density] = isa(altitude,'geopotential')

    Inputs:
    altitude -> [float] altitude in [km]

    Parameters:
    altitude_type -> [optional, string, default='geometric'] type of altitude. It can either be 'geopotential' or 'geometric'.
    Relationship between the 'geopotential' and 'geometric' altitude can be found in http://www.pdas.com/hydro.pdf

    Outputs:
    temperature -> [float] temperature in [K] at the local altitude
    pressure -> [float] pressure in [Pa] at the local altitude
    density -> [float] density in [kg/m^3] at the local altitude

    Note: the geopotential altitude should be in [610,84852] m and the geometric altitude should be in [-611,86000] m

    Reference: U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, D.C. 
               Public Domain Aeronautical Software(http://www.pdas.com/atmos.html) 
               https://gist.github.com/buzzerrookie/5b6438c603eabf13d07e
               https://ww2.mathworks.cn/help/aerotbx/ug/atmosisa.html
    '''
    altitude = altitude*1e3 # convert [km] to [m]
    # the lower atmosphere below 86km is separated into seven layers 
    geopotential_alt = [-610,11000, 20000, 32000, 47000, 51000,71000,84852] # Geopotential altitudes above MSL in [m]
    geometric_alt = [-611,11019, 20063, 32162, 47350, 51413, 71802,86000] # Geometric altitudes above MSL in [m]
    lr = np.array([-6.5, 0, 1, 2.8, 0, -2.8, -2])/1e3 # Lapse rate in [K/m]
    
    t0,p0,h0 = t_msl,p_msl,h_msl

    if altitude_type == 'geopotential':
        if altitude < geopotential_alt[0] or altitude > geopotential_alt[-1]:
            raise Exception("geopotential altitude should be in [-0.610, 84.852] km")
    elif altitude_type == 'geometric':   
        if altitude < geometric_alt[0] or altitude > geometric_alt[-1]:
            raise Exception("geometric altitude should be in [-0.611, 86.0] km") 
        # convert geometric altitude to geopotential altitude   
        altitude = altitude*R_e/(altitude+R_e)
    else:
        raise Exception("The type of altitude should either be 'geopotential' or 'geometric'")       

    for i in range(len(lr)):
        if altitude <= geopotential_alt[i+1]:
            temperature, pressure = lapse_tp(t0, p0, lr[i], h0, altitude)
            break
        else:
            # if altitudes are greater than the first several layers, then it has to integeate these layers first.
            t0, p0 = lapse_tp(t0, p0, lr[i], h0, geopotential_alt[i+1])
            h0 = geopotential_alt[i+1]

    density = pressure / (R * temperature)
    return {'temperature[K]':temperature,'pressure[Pa]':pressure,'density[kg/m^3]':density}