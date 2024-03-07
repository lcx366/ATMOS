import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import lpmv
import pkg_resources

def nrlmsis00_data():
    '''
    Read the data block from nrlmsis00_data.npz
    ''' 
    data_path = pkg_resources.resource_filename('pyatmos', 'data/')
    data = np.load(data_path+'nrlmsis00_data.npz')
    pt,pd,ps,pdl = data['pt'],data['pd'],data['ps'],data['pdl']
    ptm,pdm,ptl,pma = data['ptm'],data['pdm'],data['ptl'],data['pma']
    sam,pavgm = data['sam'],data['pavgm']
    return pt,pd,ps,pdl,ptm,pdm,ptl,pma,sam,pavgm

def tselec(switches):
    flags = {'sw':np.zeros(23),'swc':np.zeros(23)}
    for i in range(23):
        if i != 8:
            if switches[i] == 1:
                flags['sw'][i] = 1
            else:
                flags['sw'][i] = 0
            if switches[i] > 0:
                flags['swc'][i] = 1
            else:
                flags['swc'][i] = 0
        else:
            flags['sw'][i] = switches[i]
            flags['swc'][i] = switches[i]
    return flags  

def glatf(lat):
    c2 = np.cos(2*np.deg2rad(lat))
    gv = 980.616*(1 - 0.0026373*c2)
    reff = 2*gv/(3.085462E-6 + 2.27E-9*c2)*1E-5
    return gv,reff

def ccor(alt,r,h1,zh):
    e = (alt - zh)/h1
    if e > 70:
        return 1
    elif e < -70:
        return np.exp(r)
    else:
        return np.exp(r/(1 + np.exp(e)))

def ccor2(alt,r,h1,zh,h2):     
    e1 = (alt - zh)/h1
    e2 = (alt - zh)/h2
    if e1 > 70 or e2 > 70:
        return 1
    if e1 < -70 and e2 < -70:
        return np.exp(r)
    ex1,ex2 = np.exp([e1,e2])
    ccor2v = r/(1 + 0.5*(ex1 + ex2))
    return np.exp(ccor2v)

def scalh(alt,xm,temp,gsurf,re):
    rgas = 831.4
    g = rgas*temp/(gsurf/(1 + alt/re)**2*xm)
    return g

def dnet(dd,dm,zhm,xmm,xm):
    '''
    Turbopause correction for msis models
    '''      
    a  = zhm/(xmm - xm)
    if not (dm > 0 and dd > 0):
        print('dnet log error {0:.1f} {1:.1f} {2:.1f}'.format(dm,dd,xm))
        if dd == 0 and dm == 0: dd = 1
        if dm == 0: return dd
        if dd == 0: return dm
    ylog = a*np.log(dm/dd)
    if ylog < -10: return dd
    if ylog > 10: return dm
    a = dd*(1 + np.exp(ylog))**(1/a)
    return a

def zeta(zz,zl,re): 
    return (zz - zl)*(re + zl)/(re + zz)

def densm(alt, d0, xm, tz, zn3, tn3, tgn3, zn2, tn2, tgn2,gsurf,re):
    rgas = 831.4
    densm_tmp = d0
    tz_tmp = tz
    mn3,mn2 = len(zn3),len(zn2)

    if alt > zn2[0]:
        if xm == 0: 
            densm_tmp = tz
            return densm_tmp,tz_tmp       
        else:
            densm_tmp = d0
            return densm_tmp,tz_tmp 

    # stratosphere/mesosphere temperature
    if alt > zn2[mn2-1]:
        z = alt
    else:
        z = zn2[mn2-1]
    mn = mn2
    xs,ys = [np.zeros(mn) for i in range(2)]
    z1,z2 = zn2[0],zn2[mn-1]
    t1,t2=tn2[0],tn2[mn-1]
    zg,zgdif = zeta(z,z1,re),zeta(z2,z1,re)
    
    # set up spline nodes
    for k in range(mn):
        xs[k] = zeta(zn2[k],z1,re)/zgdif
        ys[k] = 1/tn2[k]
    yd1 = -tgn2[0]/t1**2*zgdif
    yd2 = -tgn2[1]/t2**2*zgdif*((re + z2)/(re + z1))**2

    # calculate spline coefficients
    cs = CubicSpline(xs,ys,bc_type=((1,yd1),(1,yd2))) 
    x = zg/zgdif
    y = cs(x)

    # temperature at altitude
    tz_tmp = 1/y
    if xm != 0:
        # calaculate stratosphere/mesospehere density
        glb = gsurf/(1 + z1/re)**2
        gamm = xm*glb*zgdif/rgas
    
        # Integrate temperature profile
        yi = cs.integrate(xs[0],x)
        expl = gamm*yi
        if expl > 50:
            expl = 50
        # Density at altitude
        densm_tmp = densm_tmp*(t1/tz_tmp)*np.exp(-expl)
    if alt > zn3[0]:
        if xm == 0:
            densm_tmp = tz_tmp
            return densm_tmp,tz_tmp      
        else:
            return densm_tmp,tz_tmp 

    # troposhere/stratosphere temperature
    z = alt
    mn = mn3
    xs,ys = [np.zeros(mn) for i in range(2)]
    z1,z2 = zn3[0],zn3[mn-1]
    t1,t2 = tn3[0],tn3[mn-1]
    zg,zgdif = zeta(z,z1,re),zeta(z2,z1,re)

    # set up spline nodes
    for k in range(mn):
        xs[k] = zeta(zn3[k],z1,re)/zgdif
        ys[k] = 1/tn3[k]
    yd1 = -tgn3[0]/t1**2*zgdif
    yd2 = -tgn3[1]/t2**2*zgdif*((re+z2)/(re+z1))**2

    # calculate spline coefficients
    cs = CubicSpline(xs,ys,bc_type=((1,yd1),(1,yd2))) 
    x = zg/zgdif
    y = cs(x)

    # temperature at altitude
    tz_tmp = 1/y
    if xm != 0:
        # calaculate tropospheric / stratosphere density
        glb = gsurf/(1 + z1/re)**2
        gamm = xm*glb*zgdif/rgas
    
        # Integrate temperature profile
        yi = cs.integrate(xs[0],x)
        expl = gamm*yi
        if expl > 50: expl = 50
        # Density at altitude
        densm_tmp = densm_tmp*(t1/tz_tmp)*np.exp(-expl)

    if xm == 0:
        densm_tmp = tz_tmp
        return densm_tmp,tz_tmp
    else:
        return densm_tmp,tz_tmp

def densu (alt,dlb,tinf,tlb,xm,alpha,tz,zlb,s2,zn1,tn1,tgn1,gsurf,re):
    '''
    Calculate Temperature and Density Profiles for MSIS models
    '''
    rgas = 831.4
    densu_tmp = 1
    mn1 = len(zn1)
    # joining altitudes of Bates and spline
    za = zn1[0]
    if alt > za:
        z = alt
    else:
        z = za
    # geopotential altitude difference from ZLB
    zg2 = zeta(z,zlb,re)

    # Bates temperature
    tt = tinf - (tinf - tlb)*np.exp(-s2*zg2)
    ta = tz = tt
    densu_tmp = tz_tmp = tz

    if alt < za:
        # calculate temperature below ZA
        # temperature gradient at ZA from Bates profile
        dta = (tinf - ta)*s2*((re + zlb)/(re + za))**2
        tgn1[0],tn1[0] = dta,ta
        if alt > zn1[mn1-1]:
            z = alt
        else:
            z = zn1[mn1-1]
        mn = mn1
        xs,ys = [np.zeros(mn) for i in range(2)]
        z1,z2 = zn1[0],zn1[mn-1]
        t1,t2 = tn1[0],tn1[mn-1]
        # geopotental difference from z1
        zg,zgdif = zeta(z,z1,re),zeta(z2,z1,re)
        # set up spline nodes
        for k in range(mn):
            xs[k] = zeta(zn1[k],z1,re)/zgdif
            ys[k] = 1/tn1[k]
        # end node derivatives
        yd1 = -tgn1[0]/t1**2*zgdif
        yd2 = -tgn1[1]/t2**2*zgdif*((re + z2)/(re + z1))**2
        # calculate spline coefficients
        cs = CubicSpline(xs,ys,bc_type=((1,yd1),(1,yd2))) 
        x = zg/zgdif
        y = cs(x)
        # temperature at altitude
        tz_tmp = 1/y
        densu_tmp = tz_tmp
    if xm == 0: return densu_tmp,tz_tmp
    
    # calculate density above za
    glb = gsurf/(1 + zlb/re)**2
    gamma = xm*glb/(s2*rgas*tinf)
    expl = np.exp(-s2*gamma*zg2)
    if expl > 50: expl = 50
    if tt <= 0: expl = 50   

    # density at altitude
    densa = dlb*(tlb/tt)**(1 + alpha + gamma)*expl
    densu_tmp = densa
    if alt >= za: return densu_tmp,tz_tmp
    
    # calculate density below za
    glb = gsurf/(1 + z1/re)**2
    gamm = xm*glb*zgdif/rgas

    # integrate spline temperatures
    yi = cs.integrate(xs[0],x)
    expl = gamm*yi
    if expl > 50: expl = 50
    if tz_tmp <= 0: expl = 50

    # density at altitude
    densu_tmp = densu_tmp*(t1/tz_tmp)**(1 + alpha)*np.exp(-expl)
    return densu_tmp,tz_tmp

# =============== 3hr magnetic activity functions =================== #

# Eq. A24d
def g0(a,p):
    return (a - 4 + (p[25] - 1)*(a - 4 + (np.exp(-np.abs(p[24])*(a - 4)) - 1) / np.abs(p[24])))

# Eq. A24c
def sumex(ex):
    return (1 + (1 - ex**19)/(1 - ex)*ex**0.5)

# Eq. A24a
def sg0(ex,p,ap):
    # call sumex, g0
    return (g0(ap[1],p) + g0(ap[2],p)*ex + g0(ap[3],p)*ex**2 + \
                g0(ap[4],p)*ex**3 + (g0(ap[5],p)*ex**4 + \
                g0(ap[6],p)*ex**12)*(1-ex**8)/(1-ex))/sumex(ex)

# =============== 3hr magnetic activity functions =================== #
def PLegendreA(lmax, x):
    """
    Calculates all unnormalized associated Legendre polynomials up to degree lmax for a given x value.

    Inputs:
        lmax -> [int] The maximum degree of the polynomials to compute.
        x -> [float] The value at which to evaluate the polynomials.
    Outputs:
        res -> [array-like] An array of shape (lmax + 1) * (lmax + 2) // 2 containing the values of the unnormalized associated Legendre polynomials P_l^m(x) for l = 0, 1, ..., lmax and m = 0, 1, ..., l.
    Note: The Condon-Shortley phase is excluded.    
    """
    res = np.zeros((lmax + 1) * (lmax + 2) // 2)

    # Compute P_l^m(x) for all l and m
    k = 0
    for l in range(lmax + 1):
        for m in range(l + 1):
            res_k = lpmv(m, l, x)
            if m%2: res_k *= -1
            res[k] = res_k
            k += 1
    return res

def lengendre(g_lat,lmax = 8):
    x = np.sin(np.deg2rad(g_lat))
    PLegendreA_x = PLegendreA(lmax,x)
    return PLegendreA_x

def globe7(p,inputp,flags):
    '''
    Calculate G(L) function 
    '''
    t = np.zeros(15)
    sr = 7.2722E-5
    dr = 1.72142E-2
    hr = 0.2618
    
    apdf = 0
    apt = np.zeros(4)
    tloc = inputp['lst']

    if not (flags['sw'][6]==0 and flags['sw'][7]==0 and flags['sw'][13]==0):
        stloc,ctloc = np.sin(hr*tloc),np.cos(hr*tloc)
        s2tloc,c2tloc = np.sin(2*hr*tloc),np.cos(2*hr*tloc)
        s3tloc,c3tloc = np.sin(3*hr*tloc),np.cos(3*hr*tloc)
    cd32 = np.cos(dr*(inputp['doy'] - p[31]))
    cd18 = np.cos(2*dr*(inputp['doy'] - p[17]))
    cd14 = np.cos(dr*(inputp['doy'] - p[13]))
    cd39 = np.cos(2*dr*(inputp['doy'] - p[38]))

    # F10.7 effect 
    df = inputp['f107'] - inputp['f107A']
    dfa = inputp['f107A'] - 150
    t[0] =  p[19]*df*(1 + p[59]*dfa) + p[20]*df**2 + p[21]*dfa + p[29]*dfa**2
    f1 = 1 + (p[47]*dfa + p[19]*df + p[20]*df**2)*flags['swc'][0]
    f2 = 1 + (p[49]*dfa + p[19]*df + p[20]*df**2)*flags['swc'][0]
    
    plg = lengendre(inputp['g_lat'])

    #  time independent 
    t[1] = p[1]*plg[3] + p[2]*plg[10] + p[22]*plg[21] + p[14]*plg[3]*dfa*flags['swc'][0] + p[26]*plg[1]
    
    # symmetrical annual
    t[2] = p[18]*cd32

    # symmetrical semiannual
    t[3] = (p[15] + p[16]*plg[3])*cd18

    # asymmetrical annual
    t[4] = f1*(p[9]*plg[1] + p[10]*plg[6])*cd14

    # asymmetrical semiannual 
    t[5] = p[37]*plg[1]*cd39
    
    # diurnal 
    if flags['sw'][6]:
        t71 = p[11]*plg[4]*cd14*flags['swc'][4]
        t72 = p[12]*plg[4]*cd14*flags['swc'][4]
        t[6] = f2*((p[3]*plg[2] + p[4]*plg[7] + p[27]*plg[16] + t71) * ctloc + (p[6]*plg[2] + p[7]*plg[7] + p[28]*plg[16] + t72)*stloc)
    
    # semiannual 
    if flags['sw'][7]:
        t81 = (p[23]*plg[8] + p[35]*plg[17])*cd14*flags['swc'][4]
        t82 = (p[33]*plg[8] + p[36]*plg[17])*cd14*flags['swc'][4]
        t[7] = f2*((p[5]*plg[5] + p[41]*plg[12] + t81)*c2tloc +(p[8]*plg[5] + p[42]*plg[12] + t82)*s2tloc)

    # terdiurnal
    if flags['sw'][13]:
        t[13] = f2*((p[39]*plg[9] + (p[93]*plg[13] + p[46]*plg[24])*cd14*flags['swc'][4])*s3tloc + (p[40]*plg[9]+(p[94]*plg[13] + p[48]*plg[24])*cd14*flags['swc'][4])*c3tloc)
    
    # magnetic activity based on daily ap 
    if flags['sw'][8] == -1:
        ap = inputp['ap_a']
        if p[51]!= 0:
            exp1 = np.exp(-10800*np.abs(p[51])/(1 + p[138]*(45 - np.abs(inputp['g_lat']))))
            if exp1 > 0.99999: exp1 = 0.99999
            if p[24] < 1E-4: p[24] = 1E-4
            apt[0] = sg0(exp1,p,ap)
            # apt[1] = sg2(exp1,p,ap)
            # apt[2] = sg0(exp2,p,ap)
            # apt[3] = sg2(exp2,p,ap)

            if flags['sw'][8]:
                t[8] = apt[0]*(p[50] + p[96]*plg[3] + p[54]*plg[10] + \
                       (p[125]*plg[1] + p[126]*plg[6] + p[127]*plg[15])*cd14*flags['swc'][4] + \
                       (p[128]*plg[2] + p[129]*plg[7] + p[130]*plg[16])*flags['swc'][6]*np.cos(hr*(tloc - p[131])))
    else:
        apd = inputp['ap'] - 4
        p44 = p[43]
        p45 = p[44]
        if p44 < 0: p44 = 1E-5
        apdf = apd + (p45 - 1)*(apd + (np.exp(-p44*apd) - 1)/p44)
        if flags['sw'][8]:
            t[8]=apdf*(p[32] + p[45]*plg[3] + p[34]*plg[10] + \
     (p[100]*plg[1] + p[101]*plg[6] + p[102]*plg[15])*cd14*flags['swc'][4] +
     (p[121]*plg[2] + p[122]*plg[7] + p[123]*plg[16])*flags['swc'][6]*np.cos(hr*(tloc - p[124])))

    if flags['sw'][9] and inputp['g_lon'] > -1000:
        # longitudinal
        if flags['sw'][10]:
            t[10] = (1 + p[80]*dfa*flags['swc'][0])*((p[64]*plg[4] + p[65]*plg[11] + p[66]*plg[22]\
                    + p[103]*plg[2] + p[104]*plg[7] + p[105]*plg[16]\
                    + flags['swc'][4]*(p[109]*plg[2] + p[110]*plg[7] + p[111]*plg[16])*cd14)*np.cos(np.deg2rad(inputp['g_lon'])) \
                    +(p[90]*plg[4]+p[91]*plg[11]+p[92]*plg[22] + p[106]*plg[2]+p[107]*plg[7]+p[108]*plg[16]\
                    + flags['swc'][4]*(p[112]*plg[2] + p[113]*plg[7] + p[114]*plg[16])*cd14)*np.sin(np.deg2rad(inputp['g_lon'])))

        # ut and mixed ut, longitude 
        if flags['sw'][11]:
            t[11]=(1 + p[95]*plg[1])*(1 + p[81]*dfa*flags['swc'][0])*\
            (1 + p[119]*plg[1]*flags['swc'][4]*cd14)*\
            ((p[68]*plg[1] + p[69]*plg[6] + p[70]*plg[15])*np.cos(sr*(inputp['sec'] - p[71])))
            t[11] += flags['swc'][10]*(p[76]*plg[8] + p[77]*plg[17] + p[78]*plg[30])*\
            np.cos(sr*(inputp['sec'] - p[79]) + 2*np.deg2rad(inputp['g_lon']))*(1 + p[137]*dfa*flags['swc'][0])
            
        # ut, longitude magnetic activity 
        if flags['sw'][10]:
            if flags['sw'][8] == -1:
                if p[51]:
                    t[12] = apt[0]*flags['swc'][10]*(1 + p[132]*plg[1])*\
                    ((p[52]*plg[4] + p[98]*plg[11] + p[67]*plg[22])* np.cos(np.deg2rad(inputp['g_lon'] - p[97])))\
                    + apt[0]*flags['swc'][10]*flags['swc'][4]*(p[133]*plg[2] + p[134]*plg[7] + p[135]*plg[16])*\
                    cd14*np.cos(np.deg2rad(inputp['g_lon'] - p[136])) + apt[0]*flags['swc'][11]* \
                    (p[55]*plg[1] + p[56]*plg[6] + p[57]*plg[15])*np.cos(sr*(inputp['sec'] - p[58]))
            else:
                t[12] = apdf*flags['swc'][10]*(1 + p[120]*plg[1])*((p[60]*plg[4] + p[61]*plg[11] + p[62]*plg[22])*\
                    np.cos(np.deg2rad(inputp['g_lon']-p[63])))+apdf*flags['swc'][10]*flags['swc'][4]* \
                    (p[115]*plg[2] + p[116]*plg[7] + p[117]*plg[16])* \
                    cd14*np.cos(np.deg2rad(inputp['g_lon'] - p[118])) \
                    + apdf*flags['swc'][11]*(p[83]*plg[1] + p[84]*plg[6] + p[85]*plg[15])* np.cos(sr*(inputp['sec'] - p[75]))

    # parms not used: 82, 89, 99, 139-149 
    tinf = p[30]
    for i in range(14):
        tinf = tinf + np.abs(flags['sw'][i])*t[i]    
    return tinf,[dfa,plg,ctloc,stloc,c2tloc,s2tloc,s3tloc,c3tloc,apdf,apt]   

def glob7s(p,inputp,flags,varli):
    pset = 2
    t = np.zeros(14)
    dr = 1.72142E-2
    [dfa,plg,ctloc,stloc,c2tloc,s2tloc,s3tloc,c3tloc,apdf,apt] = varli
    
    # confirm parameter set
    if p[99] == 0: p[99] = pset
    if p[99] != pset:
        print("Wrong parameter set for glob7s")
        return -1

    for j in range(14):
        t[j] = 0
        cd32 = np.cos(dr*(inputp['doy'] - p[31]))
        cd18 = np.cos(2*dr*(inputp['doy'] - p[17]))
        cd14 = np.cos(dr*(inputp['doy'] - p[13]))
        cd39 = np.cos(2*dr*(inputp['doy'] - p[38]))

    # F10.7 
    t[0] = p[21]*dfa

    # time independent 
    t[1] = p[1]*plg[3] + p[2]*plg[10] + p[22]*plg[21] + p[26]*plg[1] + p[14]*plg[6] + p[59]*plg[15]

    # symmetrical annual
    t[2] = (p[18] + p[47]*plg[3] + p[29]*plg[10])*cd32

    # symmetrical semiannual
    t[3] = (p[15] + p[16]*plg[3] + p[30]*plg[10])*cd18

    # asymmetrical annual
    t[4] = (p[9]*plg[1] + p[10]*plg[6] + p[20]*plg[15])*cd14

    # asymmetrical semiannual
    t[5] = p[37]*plg[1]*cd39;

    # diurnal
    if flags['sw'][6]:
        t71 = p[11]*plg[4]*cd14*flags['swc'][4]
        t72 = p[12]*plg[4]*cd14*flags['swc'][4]
        t[6] = ((p[3]*plg[2] + p[4]*plg[7] + t71)*ctloc + (p[6]*plg[2] + p[7]*plg[7] + t72)*stloc) 

    # semidiurnal
    if flags['sw'][7]:
        t81 = (p[23]*plg[8] + p[35]*plg[17])*cd14*flags['swc'][4]
        t82 = (p[33]*plg[8] + p[36]*plg[17])*cd14*flags['swc'][4]
        t[7] = ((p[5]*plg[5] + p[41]*plg[12] + t81)*c2tloc + (p[8]*plg[5] + p[42]*plg[12] + t82)*s2tloc)

    # terdiurnal
    if flags['sw'][13]:
        t[13] = p[39]*plg[9]*s3tloc + p[40]*plg[9]*c3tloc

    # magnetic activity
    if flags['sw'][8]:
        if flags['sw'][8]==1:
            t[8] = apdf * (p[32] + p[45]*plg[3]*flags['swc'][1])
        if flags['sw'][8]==-1:
            t[8]=(p[50]*apt[0] + p[96]*plg[3]*apt[0]*flags['swc'][1])

    # longitudinal
    if not (flags['sw'][9]==0 or flags['sw'][10]==0 or inputp['g_lon']<=-1000):
        t[10] = (1 + plg[1]*(p[80]*flags['swc'][4]*np.cos(dr*(inputp['doy'] - p[81]))\
            + p[85]*flags['swc'][5]*np.cos(2*dr*(inputp['doy'] - p[86])))\
            + p[83]*flags['swc'][2]*np.cos(dr*(inputp['doy'] - p[84]))\
            + p[87]*flags['swc'][3]*np.cos(2*dr*(inputp['doy'] - p[88])))\
            *((p[64]*plg[4] + p[65]*plg[11] + p[66]*plg[22]\
            + p[74]*plg[2] + p[75]*plg[7] + p[76]*plg[16])*np.cos(np.deg2rad(inputp['g_lon']))\
            + (p[90]*plg[4] + p[91]*plg[11] + p[92]*plg[22]\
            + p[77]*plg[2] + p[78]*plg[7] + p[79]*plg[16])*np.sin(np.deg2rad(inputp['g_lon'])))
    
    tt = 0
    for i in range(14):
        tt += np.abs(flags['sw'][i])*t[i]
    return tt
              
def gtd7(inputp,switches):
    tz = 0
    zn3 = np.array([32.5,20.0,15.0,10.0,0.0])
    zn2 = np.array([72.5,55.0,45.0,32.5])
    zmix= 62.5
    
    output = {'d':{'He':0,'O':0,'N2':0,'O2':0,'AR':0,'RHO':0,'H':0,'N':0,'ANM O':0},\
              't':{'TINF':0,'TG':0}}
    
    flags = tselec(switches)
    
    # Latitude variation of gravity (none for sw[1]=0) 
    xlat = inputp['g_lat']
    if flags['sw'][1]==0: xlat = 45
    gsurf,re = glatf(xlat)
    pt,pd,ps,pdl,ptm,pdm,ptl,pma,sam,pavgm = nrlmsis00_data()
    xmm = pdm[2,4]
    
    # thermosphere/mesosphere (above zn2[0])
    if inputp['alt'] > zn2[0]:
        altt = inputp['alt']
    else:
        altt = zn2[0]

    tmp = inputp['alt']
    inputp['alt'] = altt
    soutput,dm28,[meso_tn1,meso_tn2,meso_tn3,meso_tgn1,meso_tgn2,meso_tgn3],[dfa,plg,ctloc,stloc,c2tloc,s2tloc,s3tloc,c3tloc,apdf,apt] = gts7(inputp,flags,gsurf,re)
    altt = inputp['alt']
    inputp['alt'] = tmp
    # metric adjustment 
    dm28m = dm28*1E6
    output['t']['TINF'] = soutput['t']['TINF']
    output['t']['TG'] = soutput['t']['TG']
    if inputp['alt'] >= zn2[0]:
        output['d'] = soutput['d']
        return output

    varli = [dfa,plg,ctloc,stloc,c2tloc,s2tloc,s3tloc,c3tloc,apdf,apt]

    meso_tgn2[0] = meso_tgn1[1]
    meso_tn2[0] = meso_tn1[4]
    meso_tn2[1] = pma[0,0]*pavgm[0]/(1-flags['sw'][19]*glob7s(pma[0], inputp, flags,varli))
    meso_tn2[2] = pma[1,0]*pavgm[1]/(1-flags['sw'][19]*glob7s(pma[1], inputp, flags,varli))
    meso_tn2[3] = pma[2,0]*pavgm[2]/(1-flags['sw'][19]*flags['sw'][21]*glob7s(pma[2], inputp, flags,varli))
    meso_tgn2[1] = pavgm[8]*pma[9,0]*(1+flags['sw'][19]*flags['sw'][21]*glob7s(pma[9], inputp, flags,varli))*meso_tn2[3]*meso_tn2[3]/(pma[2,0]*pavgm[2])**2
    meso_tn3[0] = meso_tn2[3]
    
    if inputp['alt'] <= zn3[0]:

        meso_tgn3[0] = meso_tgn2[1]
        meso_tn3[1] = pma[3,0]*pavgm[3]/(1-flags['sw'][21]*glob7s(pma[3], inputp, flags,varli))
        meso_tn3[2] = pma[4,0]*pavgm[4]/(1-flags['sw'][21]*glob7s(pma[4], inputp, flags,varli))
        meso_tn3[3] = pma[5,0]*pavgm[5]/(1-flags['sw'][21]*glob7s(pma[5], inputp, flags,varli))
        meso_tn3[4] = pma[6,0]*pavgm[6]/(1-flags['sw'][21]*glob7s(pma[6], inputp, flags,varli))
        meso_tgn3[1] = pma[7,0]*pavgm[7]*(1+flags['sw'][21]*glob7s(pma[7], inputp, flags,varli)) *meso_tn3[4]*meso_tn3[4]/(pma[6,0]*pavgm[6])**2

    # linear transition to full mixing below znz[0]

    dmc = 0
    if inputp['alt'] > zmix:
        dmc = 1 - (zn2[0]-inputp['alt'])/(zn2[0] - zmix)
    dz28 = soutput['d']['N2']
    
    # N2 density
    dmr = soutput['d']['N2'] / dm28m - 1
    output['d']['N2'],tz = densm(inputp['alt'],dm28m,xmm, tz, zn3, meso_tn3, meso_tgn3, zn2, meso_tn2, meso_tgn2,gsurf,re)
    output['d']['N2'] = output['d']['N2'] * (1 + dmr*dmc)

    # HE density 
    dmr = soutput['d']['He'] / (dz28 * pdm[0,1]) - 1
    output['d']['He'] = output['d']['N2'] * pdm[0,1] * (1 + dmr*dmc)

    # O density
    output['d']['O'] = 0
    output['d']['ANM O'] = 0

    # O2 density
    dmr = soutput['d']['O2'] / (dz28 * pdm[3,1]) - 1
    output['d']['O2'] = output['d']['N2'] * pdm[3,1] * (1 + dmr*dmc)

    # AR density 
    dmr = soutput['d']['AR'] / (dz28 * pdm[4,1]) - 1
    output['d']['AR'] = output['d']['N2'] * pdm[4,1] * (1 + dmr*dmc)

    # Hydrogen density
    output['d']['H'] = 0

    # Atomic nitrogen density 
    output['d']['N'] = 0

    # Total mass density 
    output['d']['RHO'] = 1.66E-24 * (4 * output['d']['He'] + 16 * output['d']['O'] + 28 * output['d']['N2']\
                                     + 32 * output['d']['O2'] + 40 * output['d']['AR'] + output['d']['H'] + 14 * output['d']['N'])

    output['d']['RHO'] = output['d']['RHO']/1000

    # temperature at altitude 
    dd,tz = densm(inputp['alt'], 1, 0, tz, zn3, meso_tn3, meso_tgn3, zn2, meso_tn2, meso_tgn2,gsurf,re)
    output['t']['TG'] = tz
    return output


def gtd7d(inputp, flags):
    output = gtd7(inputp, flags)
    output['d']['RHO'] = 1.66E-24 * (4 * output['d']['He'] + 16 * output['d']['O'] + 28 * output['d']['N2']\
                                     + 32 * output['d']['O2'] + 40 * output['d']['AR'] + output['d']['H'] + 14 * output['d']['N'] + 16 * output['d']['ANM O'])

    output['d']['RHO'] = output['d']['RHO']/1e3
    return output
              
def gts7(inputp,flags,gsurf,re):
    
    output = {'d':{'He':0,'O':0,'N2':0,'O2':0,'AR':0,'RHO':0,'H':0,'N':0,'ANM O':0},\
              't':{'TINF':0,'TG':0}}
    tz = 0
    dm28 = 0
    meso_tn1,meso_tn3 = [np.zeros(5) for i in range(2)]
    meso_tn2 = np.zeros(4)
    meso_tgn1,meso_tgn2,meso_tgn3 = [np.zeros(2) for i in range(3)]
    
    zn1 = np.array([120.0, 110.0, 100.0, 90.0, 72.5])

    dr = 1.72142E-2
    alpha = np.array([-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0])
    altl = np.array([200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0])
    pt,pd,ps,pdl,ptm,pdm,ptl,pma,sam,pavgm = nrlmsis00_data()
    za = pdl[1,15]
    zn1[0] = za
    
    # tinf variations not important below za or zn1[0]
    if inputp['alt'] > zn1[0]:
        tinf_tmp,varli = globe7(pt,inputp,flags)
        tinf = ptm[0]*pt[0] * (1+flags['sw'][15]*tinf_tmp)
    else:
        tinf = ptm[0]*pt[0]
    output['t']['TINF'] = tinf
    
    # gradient variations not important below zn1[4]
    if inputp['alt'] > zn1[4]:
        tinf_tmp,varli = globe7(ps,inputp,flags)
        grad = ptm[3]*ps[0] * (1+flags['sw'][18]*tinf_tmp)
    else:
        grad = ptm[3]*ps[0]
    tinf_tmp,varli = globe7(pd[3],inputp,flags)    
    tlb = ptm[1] * (1 + flags['sw'][16]*tinf_tmp)*pd[3,0]
    s = grad/(tinf - tlb)
    
    # Lower thermosphere temp variations not significant for density above 300 km
    if inputp['alt'] < 300:
        meso_tn1[1] = ptm[6]*ptl[0,0]/(1.0-flags['sw'][17]*glob7s(ptl[0], inputp, flags,varli))
        meso_tn1[2] = ptm[2]*ptl[1,0]/(1.0-flags['sw'][17]*glob7s(ptl[1], inputp, flags,varli))
        meso_tn1[3] = ptm[7]*ptl[2,0]/(1.0-flags['sw'][17]*glob7s(ptl[2], inputp, flags,varli))
        meso_tn1[4] = ptm[4]*ptl[3,0]/(1.0-flags['sw'][17]*flags['sw'][19]*glob7s(ptl[3], inputp, flags,varli))
        meso_tgn1[1] = ptm[8]*pma[8,0]*(1.0+flags['sw'][17]*flags['sw'][19]*glob7s(pma[8], inputp, flags,varli))*meso_tn1[4]*meso_tn1[4]/(ptm[4]*ptl[3,0])**2
    else:
        meso_tn1[1]=ptm[6]*ptl[0,0]
        meso_tn1[2]=ptm[2]*ptl[1,0]
        meso_tn1[3]=ptm[7]*ptl[2,0]
        meso_tn1[4]=ptm[4]*ptl[3,0]
        meso_tgn1[1]=ptm[8]*pma[8,0]*meso_tn1[4]*meso_tn1[4]/(ptm[4]*ptl[3,0])**2
        
    # N2 variation factor at Zlb
    tinf_tmp,varli = globe7(pd[2],inputp,flags)
    g28 = flags['sw'][20]*tinf_tmp

    # variation of turbopause height
    zhf = pdl[1,24]*(1+flags['sw'][4]*pdl[0,24]*np.sin(np.deg2rad(inputp['g_lat']))*np.cos(dr*(inputp['doy']-pt[13])))
    output['t']['TINF'] = tinf
    xmm = pdm[2,4]
    z = inputp['alt']

    # N2 density
    # Diffusive density at Zlb 
    db28 = pdm[2,0]*np.exp(g28)*pd[2,0]
    # Diffusive density at Alt 
    output['d']['N2'],output['t']['TG'] = densu(z,db28,tinf,tlb,28,alpha[2],output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
    dd = output['d']['N2']
    # Turbopause 
    zh28 = pdm[2,2]*zhf
    zhm28 = pdm[2,3]*pdl[1,5] 
    xmd = 28 - xmm
    # Mixed density at Zlb 
    b28,tz = densu(zh28,db28,tinf,tlb,xmd,(alpha[2]-1),tz,ptm[5],s, zn1,meso_tn1,meso_tgn1,gsurf,re)
    if flags['sw'][14] and z <= altl[2]:
        # Mixed density at Alt 
        dm28,tz = densu(z,b28,tinf,tlb,xmm,alpha[2],tz,ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        # Net density at Alt
        output['d']['N2'] = dnet(output['d']['N2'],dm28,zhm28,xmm,28)
    
    # HE density
    # Density variation factor at Zlb
    tinf_tmp,varli = globe7(pd[0],inputp,flags)
    g4 = flags['sw'][20]*tinf_tmp
    # Diffusive density at Zlb 
    db04 = pdm[0,0]*np.exp(g4)*pd[0,0]
    # Diffusive density at Alt 
    output['d']['He'],output['t']['TG'] = densu(z,db04,tinf,tlb, 4,alpha[0],output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
    dd = output['d']['He']
    if flags['sw'][14] and z<altl[0]:
        # Turbopause 
        zh04 = pdm[0,2]
        # Mixed density at Zlb
        b04,output['t']['TG'] = densu(zh04,db04,tinf,tlb,4-xmm,alpha[0]-1,output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        # Mixed density at Alt
        dm04,output['t']['TG'] = densu(z,b04,tinf,tlb,xmm,0,output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        zhm04 = zhm28
        # Net density at Alt
        output['d']['He'] = dnet(output['d']['He'],dm04,zhm04,xmm,4)
        # Correction to specified mixing ratio at ground 
        rl = np.log(b28*pdm[0,1]/b04)
        zc04 = pdm[0,4]*pdl[1,0]
        hc04 = pdm[0,5]*pdl[1,1]
        # Net density corrected at Alt 
        output['d']['He'] = output['d']['He']*ccor(z,rl,hc04,zc04) 
        
    # O density
    # Density variation factor at Zlb 
    tinf_tmp,varli = globe7(pd[1],inputp,flags)
    g16 = flags['sw'][20]*tinf_tmp
    # Diffusive density at Zlb 
    db16 =  pdm[1,0]*np.exp(g16)*pd[1,0]
    # Diffusive density at Alt 
    output['d']['O'],output['t']['TG'] = densu(z,db16,tinf,tlb,16,alpha[1],output['t']['TG'],ptm[5],s, zn1,meso_tn1,meso_tgn1,gsurf,re)
    dd = output['d']['O']
    if flags['sw'][14] and z <= altl[1]:
        #   Turbopause 
        zh16 = pdm[1,2]
        #  Mixed density at Zlb 
        b16,output['t']['TG'] = densu(zh16,db16,tinf,tlb,16-xmm,alpha[1]-1, output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        #  Mixed density at Alt 
        dm16,output['t']['TG'] = densu(z,b16,tinf,tlb,xmm,0,output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        zhm16 = zhm28
        # Net density at Alt 
        output['d']['O'] = dnet(output['d']['O'],dm16,zhm16,xmm,16)
        rl = pdm[1,1]*pdl[1,16]*(1+flags['sw'][0]*pdl[0,23]*(inputp['f107A']-150))
        hc16 = pdm[1,5]*pdl[1,3]
        zc16 = pdm[1,4]*pdl[1,2]
        hc216 = pdm[1,5]*pdl[1,4]
        output['d']['O'] = output['d']['O']*ccor2(z,rl,hc16,zc16,hc216)
        # Chemistry correction 
        hcc16 = pdm[1,7]*pdl[1,13]
        zcc16 = pdm[1,6]*pdl[1,12]
        rc16 = pdm[1,3]*pdl[1,14]
        # Net density corrected at Alt
        output['d']['O'] = output['d']['O']*ccor(z,rc16,hcc16,zcc16)
    
    # O2 density
    # Density variation factor at Zlb 
    tinf_tmp,varli = globe7(pd[4],inputp,flags)
    g32 = flags['sw'][20]*tinf_tmp
    # Diffusive density at Zlb 
    db32 = pdm[3,0]*np.exp(g32)*pd[4,0]
    # Diffusive density at Alt 
    output['d']['O2'],output['t']['TG'] = densu(z,db32,tinf,tlb, 32,alpha[3],output['t']['TG'],ptm[5],s, zn1,meso_tn1,meso_tgn1,gsurf,re)
    dd = output['d']['O2'];
    if flags['sw'][14]:
        if z <= altl[3]:
            # Turbopause 
            zh32 = pdm[3,2]
            # Mixed density at Zlb
            b32,output['t']['TG'] = densu(zh32,db32,tinf,tlb,32-xmm,alpha[3]-1, output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
            #  Mixed density at Alt 
            dm32,output['t']['TG'] = densu(z,b32,tinf,tlb,xmm,0,output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
            zhm32 = zhm28
            # Net density at Alt
            output['d']['O2'] = dnet(output['d']['O2'],dm32,zhm32,xmm,32)
            # Correction to specified mixing ratio at ground 
            rl = np.log(b28*pdm[3,1]/b32)
            hc32 = pdm[3,5]*pdl[1,7]
            zc32 = pdm[3,4]*pdl[1,6]
            output['d']['O2'] = output['d']['O2']*ccor(z,rl,hc32,zc32)

        # Correction for general departure from diffusive equilibrium above Zlb */
        hcc32 = pdm[3,7]*pdl[1,22]
        hcc232 = pdm[3,7]*pdl[0,22]
        zcc32 = pdm[3,6]*pdl[1,21]
        rc32 = pdm[3,3]*pdl[1,23]*(1+flags['sw'][0]*pdl[0,23]*(inputp['f107A']-150))
        # Net density corrected at Alt 
        output['d']['O2'] = output['d']['O2']*ccor2(z,rc32,hcc32,zcc32,hcc232)
    # AR density
    # Density variation factor at Zlb 
    tinf_tmp,varli = globe7(pd[5],inputp,flags)
    g40 = flags['sw'][20]*tinf_tmp
    # Diffusive density at Zlb 
    db40 = pdm[4,0]*np.exp(g40)*pd[5,0]
    # Diffusive density at Alt
    output['d']['AR'],output['t']['TG'] = densu(z,db40,tinf,tlb, 40,alpha[4],output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
    dd = output['d']['AR']
    if flags['sw'][14] and z <= altl[4]:
        # Turbopause
        zh40 = pdm[4,2]
        # Mixed density at Zlb 
        b40,output['t']['TG'] = densu(zh40,db40,tinf,tlb,40-xmm,alpha[4]-1,output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        # Mixed density at Alt
        dm40,output['t']['TG'] = densu(z,b40,tinf,tlb,xmm,0,output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        zhm40 = zhm28
        # Net density at Alt 
        output['d']['AR'] = dnet(output['d']['AR'],dm40,zhm40,xmm,40)
        # Correction to specified mixing ratio at ground 
        rl = np.log(b28*pdm[4,1]/b40)
        hc40 = pdm[4,5]*pdl[1,9]
        zc40 = pdm[4,4]*pdl[1,8]
        # Net density corrected at Alt
        output['d']['AR'] = output['d']['AR']*ccor(z,rl,hc40,zc40)
        
    # Hydrogen density 
    # Density variation factor at Zlb */
    tinf_tmp,varli = globe7(pd[6], inputp, flags)
    g1 = flags['sw'][20]*tinf_tmp
    # Diffusive density at Zlb 
    db01 = pdm[5,0]*np.exp(g1)*pd[6,0]
    # Diffusive density at Alt
    output['d']['H'],output['t']['TG']=densu(z,db01,tinf,tlb,1,alpha[6],output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
    dd = output['d']['H']
    if flags['sw'][14] and z <= altl[6]:
        # Turbopause 
        zh01 = pdm[5,2]
        # Mixed density at Zlb
        b01,output['t']['TG'] = densu(zh01,db01,tinf,tlb,1-xmm,alpha[6]-1, output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        # Mixed density at Alt 
        dm01,output['t']['TG'] = densu(z,b01,tinf,tlb,xmm,0,output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        zhm01 = zhm28
        # Net density at Alt
        output['d']['H'] = dnet(output['d']['H'],dm01,zhm01,xmm,1)
        # Correction to specified mixing ratio at ground 
        rl = np.log(b28*pdm[5,1]*np.abs(pdl[1,17])/b01)
        hc01 = pdm[5,5]*pdl[1,11]
        zc01 = pdm[5,4]*pdl[1,10]
        output['d']['H'] = output['d']['H']*ccor(z,rl,hc01,zc01)
        # Chemistry correction 
        hcc01 = pdm[5,7]*pdl[1,19]
        zcc01 = pdm[5,6]*pdl[1,18]
        rc01 = pdm[5,3]*pdl[1,20]
        # Net density corrected at Alt
        output['d']['H'] = output['d']['H']*ccor(z,rc01,hcc01,zcc01)
    
    # Atomic Nitrogen density 
    # Density variation factor at Zlb */
    tinf_tmp,varli = globe7(pd[7],inputp,flags)
    g14 = flags['sw'][20]*tinf_tmp
    # Diffusive density at Zlb 
    db14 = pdm[6,0]*np.exp(g14)*pd[7,0]
    # Diffusive density at Alt 
    output['d']['N'],output['t']['TG']=densu(z,db14,tinf,tlb,14,alpha[7],output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
    dd = output['d']['N']
    if flags['sw'][14] and z <= altl[7]:
        # Turbopause
        zh14 = pdm[6,2]
        # Mixed density at Zlb
        b14,output['t']['TG'] = densu(zh14,db14,tinf,tlb,14-xmm,alpha[7]-1, output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        # Mixed density at Alt 
        dm14,output['t']['TG'] = densu(z,b14,tinf,tlb,xmm,0,output['t']['TG'],ptm[5],s,zn1,meso_tn1,meso_tgn1,gsurf,re)
        zhm14 = zhm28
        # Net density at Alt
        output['d']['N'] = dnet(output['d']['N'],dm14,zhm14,xmm,14)
        # Correction to specified mixing ratio at ground 
        rl = np.log(b28*pdm[6,1]*np.abs(pdl[0,2])/b14)
        hc14 = pdm[6,5]*pdl[0,1]
        zc14 = pdm[6,4]*pdl[0,0]
        output['d']['N'] = output['d']['N']*ccor(z,rl,hc14,zc14)
        # Chemistry correction
        hcc14 = pdm[6,7]*pdl[0,4]
        zcc14 = pdm[6,6]*pdl[0,3]
        rc14 = pdm[6,3]*pdl[0,5]
        # Net density corrected at Alt
        output['d']['N'] = output['d']['N']*ccor(z,rc14,hcc14,zcc14)
    
    # Anomalous Oxygen density 
    tinf_tmp,varli = globe7(pd[8],inputp,flags)
    g16h = flags['sw'][20]*tinf_tmp
    db16h = pdm[7,0]*np.exp(g16h)*pd[8,0]
    tho = pdm[7,9]*pdl[0,6]
    dd,output['t']['TG'] = densu(z,db16h,tho,tho,16,alpha[8],output['t']['TG'],ptm[5],s, zn1,meso_tn1,meso_tgn1,gsurf,re)
    zsht = pdm[7,5]
    zmho = pdm[7,4]
    zsho = scalh(zmho,16,tho,gsurf,re)
    output['d']['ANM O'] = dd*np.exp(-zsht/zsho*(np.exp(-(z-zmho)/zsht)-1))

    # total mass density
    output['d']['RHO'] = 1.66E-24*(4*output['d']['He']+16*output['d']['O']+28*output['d']['N2']\
                                   +32*output['d']['O2']+40*output['d']['AR']+ output['d']['H']+14*output['d']['N'])

    # temperature 
    z = inputp['alt']
    ddum,output['t']['TG'] = densu(z,1, tinf, tlb, 0, 0, output['t']['TG'], ptm[5], s, zn1, meso_tn1, meso_tgn1,gsurf,re)

    # convert to g/cm^3 
    for key in output['d'].keys():
        output['d'][key] = output['d'][key]*1.0E6
    output['d']['RHO'] = output['d']['RHO']/1000    
    return output,dm28,[meso_tn1,meso_tn2,meso_tn3,meso_tgn1,meso_tgn2,meso_tgn3],varli