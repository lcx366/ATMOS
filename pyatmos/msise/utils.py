# ------------------------------------------------------------------- #
# ----------------------------- utilities --------------------------- #
# ------------------------------------------------------------------- #

'''
pyatmos utils

This submodule defines the following functions:

wraplon - Wrap a longitude in range of [0,360] to [-180,180]

hms2s   - Convert hour/minute/second to seconds

hms2h   - Convert hour/minute/second to hours        
'''

# ========================= convert position ======================== #

def wraplon(lon):
    if lon > 180:
        lonwrap = lon - 360
    else:
        lonwrap = lon
    return lonwrap 

# =========================== convert time ========================== #

def hms2s(h,m,s):
    return h*3.6E3 + m*60 + s

def hms2h(h,m,s):
    return h + m/60 + s/3.6E3

