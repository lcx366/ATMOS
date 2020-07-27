def wraplon(lon):
    '''
    Wrap a longitude in range of [0,360] to [-180,180]
    '''
    if lon > 180:
        lonwrap = lon - 360
    else:
        lonwrap = lon
    return lonwrap 

def hms2s(h,m,s):
    '''
    Convert hour/minute/second to seconds
    '''
    return h*3.6E3 + m*60 + s

def hms2h(h,m,s):
    '''
    Convert hour/minute/second to hours
    '''
    return h + m/60 + s/3.6E3