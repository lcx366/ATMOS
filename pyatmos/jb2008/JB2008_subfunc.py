import numpy as np
from numba import jit

from ..utils import Const

@jit(nopython=True)
def JB2008(AMJD,YRDAY,SUN,SAT,F10,F10B,S10,S10B,M10,M10B,Y10,Y10B,DSTDTC):
    '''
    Jacchia-Bowman 2008 Model Atmosphere

    This is the CIRA "Integration Form" of a Jacchia Model.
    There are no tabular values of density. Instead, the barometricequation and diffusion equation are integrated numerically using
    the Newton-Coates method to produce the density profile up to the input position.

    INPUT:
        AMJD   : Date and Time, in modified Julian Days and Fraction (MJD = JD-2400000.5)
        SUN[0] : Right Ascension of Sun (radians)
        SUN[1] : Declination of Sun (radians)
        SAT[0] : Right Ascension of Position (radians)
        SAT[1] : Geocentric Latitude of Position (radians)
        SAT[2] : Height of Position (km)
        F10    : 10.7-cm Solar Flux (1.0E-22*W/(M**2*Hz)) (Tabular time 1.0 day earlier)
        F10B   : 10.7-cm Solar Flux, ave. 81-day centered on the input time (Tabular time 1.0 day earlier)
        S10    : EUV index (26-34 nm) scaled to F10 (Tabular time 1.0 day earlier)
        S10B   : EUV 81-day ave. centered index (Tabular time 1.0 day earlier)
        M10   : MG2 index scaled to F10 (Tabular time 2.0 days earlier)
        M10B  : MG2 81-day ave. centered index (Tabular time 2.0 days earlier)
        Y10    : Solar X-Ray & Lya index scaled to F10 (Tabular time 5.0 days earlier)
        Y10B   : Solar X-Ray & Lya 81-day ave. centered index (Tabular time 5.0 days earlier)
        DSTDTC : Temperature change computed from Dst index

    OUTPUT:
        TEMP(1): Exospheric Temperature above Input Position (K)
        TEMP(2): Temperature at Input Position (K)
        RHO    : Total Mass-Desnity at Input Position (kg/m**3)

    Reference:
        Bowman, Bruce R., etc. : "A New Empirical Thermospheric
        Density Model JB2008 Using New Solar and Geomagnetic Indices",
        AIAA/AAS 2008, COSPAR CIRA 2008 Model

    Note:
        The program is translated from the fortran source code written by Bruce R Bowman (HQ AFSPC, Space Analysis Division), 2008

    '''

    # The alpha are the thermal diffusion coefficients in Eq. (6)
    TEMP = np.zeros(2)
    ALPHA = np.zeros(5)
    ALPHA[4] = -0.38

    # AL10 is DLOG(10.0)
    AL10 = Const.al10

    # The AMW are the molecular weights in order: N2, O2, O, Ar, He & H
    AMW = np.array([28.0134,31.9988,15.9994,39.9480,4.0026,1.00797])

    # AVOGAD is Avogadro's number in mks units (molecules/kmol)
    AVOGAD = Const.avogad

    PI = Const.pi
    TWOPI,FOURPI = Const.twopi,Const.fourpi
    PIOV2,PIOV4 = Const.pivo2,Const.pivo4

    # The FRAC are the assumed sea-level volume fractions in order: N2, O2, Ar, and He
    FRAC = np.array([0.7811,0.20955,9.34e-3,1.289e-5])

    # RSTAR is the universal gas-constant in mks units (joules/K/kmol)
    RSTAR = Const.rstar

    # The R# are values used to establish height step sizes in the regimes 90km to 105km, 105km to 500km and 500km upward.
    R1,R2,R3 = 0.01,0.025,0.075

    # The WT are weights for the Newton-Cotes Five-Point Quad. formula
    WT = np.array([7,32,12,32,7])*2/45

    # The CHT are coefficients for high altitude density correction
    CHT = np.array([0.22,-0.2e-2,0.115e-2,-0.211e-5])
    DEGRAD  =  Const.degrad 

    # Equation (14)
    FN = (F10B/240)**0.25
    if FN > 1: FN = 1
    FSB = F10B*FN + S10B*(1 - FN)
    TSUBC = 392.4 + 3.227*FSB + 0.298*(F10-F10B) + 2.259*(S10-S10B) + 0.312*(M10-M10B) + 0.178*(Y10-Y10B)

    # Equation (15)
    ETA = np.abs(SAT[1] - SUN[1])/2
    THETA = np.abs(SAT[1] + SUN[1])/2

    # Equation (16)
    H = SAT[0] - SUN[0]
    TAU = H - 0.64577182 + 0.10471976 * np.sin(H + 0.75049158)

    GLAT,ZHT  = SAT[1],SAT[2]
    GLST  = H + PI
    GLSTHR = (GLST/DEGRAD)*(24/360)
    if GLSTHR >= 24: GLSTHR -= 24
    if GLSTHR < 0: GLSTHR += 24

    # Equation (17)
    C = np.cos(ETA)**2.5
    S = np.sin(THETA)**2.5

    DF = S + (C - S) * np.abs(np.cos(0.5 * TAU))**3
    TSUBL = TSUBC * (1 + 0.31 * DF)

    # Compute correction to dTc for local solar time and lat correction
    DTCLST = DTSUB(F10,GLSTHR,GLAT,ZHT)

    # Compute the local exospheric temperature.
    # Add geomagnetic storm effect from input dTc value

    TEMP[0] = TSUBL + DSTDTC
    TINF = TEMP[0] + DTCLST

    # Equation (9)
    TSUBX = 444.3807 + 0.02385 * TINF - 392.8292 * np.exp(-0.0021357 * TINF)

    # Equation (11)
    GSUBX = 0.054285714 * (TSUBX - 183)

    # The TC array will be an argument in the call to XLOCAL, which evaluates Equation (10) or Equation (13)
    TC = np.zeros(4)
    TC[0],TC[1] = TSUBX,GSUBX

    # A AND GSUBX/A OF Equation (13)
    TC[2] = (TINF - TSUBX)/PIOV2
    TC[3] = GSUBX/TC[2]
    
    # Equation (5)
    Z1 = 90
    Z2 = min(SAT[2],105)
    AL = np.log(Z2/Z1)
    N = int(AL/R1) + 1
    ZR = np.exp(AL/N)
    AMBAR1 = XAMBAR(Z1)
    TLOC1 = XLOCAL(Z1,TC)
    ZEND = Z1
    SUM2 = 0
    AIN = AMBAR1 * XGRAV(Z1)/TLOC1

    for I in range(N):
        Z = ZEND
        ZEND = ZR * Z
        DZ = 0.25 * (ZEND-Z)
        SUM1 = WT[0]*AIN

        for J in range(1,5):
            Z += DZ
            AMBAR2 = XAMBAR(Z)
            TLOC2 = XLOCAL(Z,TC)
            GRAVL = XGRAV(Z)
            AIN = AMBAR2 * GRAVL/TLOC2
            SUM1 += WT[J] * AIN
        SUM2 += DZ * SUM1
    
    FACT1 = 1e3/RSTAR
    RHO = 3.46e-6 * AMBAR2 * TLOC1 * np.exp(-FACT1*SUM2) /AMBAR1 /TLOC2

    # Equation (2)
    ANM = AVOGAD * RHO
    AN = ANM/AMBAR2

    # Equation (3)
    FACT2  = ANM/28.96
    ALN = np.zeros(6)
    ALN[0] = np.log(FRAC[0]*FACT2)
    ALN[3] = np.log(FRAC[2]*FACT2)
    ALN[4] = np.log(FRAC[3]*FACT2)

    # Equation (4)
    ALN[1] = np.log(FACT2 * (1 + FRAC[1]) - AN)
    ALN[2] = np.log(2 * (AN - FACT2))
 
    if SAT[2] <= 105:
        TEMP[1] = TLOC2
        # Put in negligible hydrogen for use in DO-LOOP 13
        ALN[5] = ALN[4] - 25
    else:    
        # Equation (6)
        Z3 = min(SAT[2],500)
        AL = np.log(Z3/Z)
        N = int(AL/R2) + 1
        ZR = np.exp(AL/N)
        SUM2 = 0
        AIN = GRAVL/TLOC2

        for I in range(N):
            Z = ZEND
            ZEND = ZR * Z
            DZ = 0.25 * (ZEND - Z)
            SUM1 = WT[0] * AIN
            for J in range(1,5):
                Z += DZ
                TLOC3 = XLOCAL(Z,TC)
                GRAVL = XGRAV(Z)
                AIN = GRAVL/TLOC3
                SUM1 += WT[J] * AIN
            SUM2 += DZ * SUM1 

        Z4 = max(SAT[2],500)
        AL = np.log(Z4/Z)
        R = R2
        if SAT[2] > 500: R = R3
        N = int(AL/R) + 1
        ZR = np.exp(AL/N)
        SUM3 = 0
        for I in range(N):
            Z = ZEND
            ZEND = ZR * Z
            DZ = 0.25 * (ZEND - Z)
            SUM1 = WT[0] * AIN
            for J in range(1,5):
                Z += DZ
                TLOC4 = XLOCAL(Z,TC)
                GRAVL = XGRAV(Z)
                AIN = GRAVL/TLOC4
                SUM1 += WT[J] * AIN
            SUM3 += DZ * SUM1
        
        if SAT[2] <= 500:
            T500 = TLOC4
            TEMP[1] = TLOC3
            ALTR = np.log(TLOC3/TLOC2)
            FACT2 = FACT1 * SUM2
            HSIGN = 1
        else:
            T500 = TLOC3
            TEMP[1] = TLOC4
            ALTR = np.log(TLOC4/TLOC2)
            FACT2 = FACT1 * (SUM2 + SUM3)
            HSIGN = -1
        
        ALN[:-1] -= (1 + ALPHA) * ALTR + FACT2 * AMW[:-1]  
        # Equation (7) - Note that in CIRA72, AL10T5 = DLOG10(T500)
        AL10T5 = np.log10(TINF)
        ALNH5 = (5.5 * AL10T5 - 39.4) * AL10T5 + 73.13
        ALN[5] = AL10 * (ALNH5 + 6) + HSIGN * (np.log(TLOC4/TLOC3) + FACT1 * SUM3 * AMW[5])

    # Equation (24)  - J70 Seasonal-Latitudinal Variation
    TRASH = (AMJD - 36204) / 365.2422
    CAPPHI = TRASH%1
    DLRSL = 0.02 * (SAT[2] - 90) * np.exp(-0.045 * (SAT[2] - 90)) * np.sign(SAT[1]) * np.sin(TWOPI * CAPPHI+ 1.72) * np.sin(SAT[1])**2

    # Equation (23) - Computes the semiannual variation
    DLRSA = 0
    if Z < 2e3:
        # Use new semiannual model
        FZZ,GTZ,DLRSA = SEMIAN08(YRDAY,ZHT,F10B,S10B,M10B)
        if FZZ < 0: DLRSA = 0  
    
    # Sum the delta-log-rhos and apply to the number densities.
    # In CIRA72 the following equation contains an actual sum, namely DLR = AL10 * (DLRGM + DLRSA + DLRSL)
    # However, for Jacchia 70, there is no DLRGM or DLRSA.

    DLR = AL10 * (DLRSL + DLRSA)
    ALN += DLR
    # Compute mass-density and mean-molecular-weight and convert number density logs from natural to common.
    AN = np.exp(ALN)   
    SUMNM = (AN*AMW).sum()
    AL10N = ALN/AL10
    RHO = SUMNM/AVOGAD
    
    # Compute the high altitude exospheric density correction factor
    FEX = 1

    if ZHT >= 1e3 and ZHT < 1.5e3:
        ZETA = (ZHT - 1e3) * 0.002
        ZETA2 =  ZETA**2
        ZETA3 =  ZETA**3
        F15C = CHT[0] + CHT[1]*F10B + (CHT[2] + CHT[3]*F10B)*1.5e3
        F15C_ZETA = (CHT[2] + CHT[3]*F10B) * 500
        FEX2 = 3 * F15C - F15C_ZETA - 3
        FEX3 = F15C_ZETA - 2 * F15C + 2
        FEX = 1 + FEX2 * ZETA2 + FEX3 * ZETA3

    if ZHT >= 1.5e3: FEX = CHT[0] + CHT[1]*F10B + CHT[2]*ZHT + CHT[3]*F10B*ZHT

    # Apply the exospheric density correction factor.
    RHO *= FEX

    return TEMP,RHO

@jit(nopython=True)
def XAMBAR(Z):
    '''
    Evaluates Equation (1)
    '''
    C = np.array([28.15204,-8.5586e-2,1.2840e-4,-1.0056e-5,-1.0210e-5,1.5044e-6,9.9826e-8])
    DZ = Z - 100
    AMB = C[6]
    for i in range(5,-1,-1): AMB = DZ * AMB + C[i]
    return AMB

@jit(nopython=True)
def XGRAV(Z):
    '''
    Evaluates Equation (8)
    '''
    return 9.80665/(1 + Z/6356.766)**2

@jit(nopython=True)
def XLOCAL(Z,TC):
    '''
    Evaluates Equation (10) or Equation (13), depending on Z
    '''
    DZ = Z - 125
    if DZ > 0:
      XLOCAL = TC[0] + TC[2] * np.arctan(TC[3]*DZ*(1 + 4.5e-6*DZ**2.5))
    else:      
      XLOCAL = ((-9.8204695e-6 * DZ - 7.3039742e-4) * DZ**2 + 1)* DZ * TC[1] + TC[0]
    return XLOCAL

@jit(nopython=True)
def DTSUB (F10,XLST,XLAT,ZHT):
    '''
    COMPUTE dTc correction for Jacchia-Bowman model

    Calling Args:

    F10       = (I)   F10 FLUX
    XLST      = (I)   LOCAL SOLAR TIME (HOURS 0-23.999)
    XLAT      = (I)   XLAT = SAT LAT (RAD)
    ZHT       = (I)   ZHT = HEIGHT (KM)
    DTC       = (O)   dTc correction
    '''
    B = np.array([-4.57512297, -5.12114909, -69.3003609,\
                   203.716701,  703.316291, -1943.49234,\
                   1106.51308, -174.378996,  1885.94601,\
                   -7093.71517,  9224.54523, -3845.08073,\
                   -6.45841789,  40.9703319, -482.006560,\
                   1818.70931, -2373.89204,  996.703815,36.1416936])
    C = np.array([-15.5986211, -5.12114909, -69.3003609,\
                  203.716701,  703.316291, -1943.49234,\
                  1106.51308, -220.835117,  1432.56989,\
                  -3184.81844, 3289.81513, -1353.32119,\
                  19.9956489, -12.7093998,  21.2825156,\
                  -2.75555432, 11.0234982,  148.881951,\
                  -751.640284, 637.876542,  12.7093998,\
                  -21.2825156, 2.75555432])
    DTC = 0
    tx = XLST/24
    ycs = np.cos(XLAT)
    F = (F10 - 100)/100

    # calculates dTc
    if ZHT >= 120 and ZHT <= 200:
        H = (ZHT - 200)/50
        DTC200 = C[16] + C[17]*tx*ycs + C[18]*tx**2*ycs + C[19]*tx**3*ycs + C[20]*F*ycs + C[21]*tx*F*ycs + C[22]*tx**2*F*ycs

        sum_ = C[0] + B[1]*F + C[2]*tx*F + C[3]*tx**2*F + C[4]*tx**3*F + C[5]*tx**4*F + C[6]*tx**5*F +\
               C[7]*tx*ycs + C[8]*tx**2*ycs + C[9]*tx**3*ycs + C[10]*tx**4*ycs + C[11]*tx**5*ycs + C[12]*ycs+\
               C[13]*F*ycs + C[14]*tx*F*ycs  + C[15]*tx**2*F*ycs

        DTC200DZ = sum_
        CC = 3*DTC200 - DTC200DZ
        DD = DTC200 - CC
        ZP = (ZHT-120)/80
        DTC = CC*ZP*ZP + DD*ZP*ZP*ZP

    if ZHT > 200 and ZHT <= 240:
        H = (ZHT - 200)/50
        sum_ = C[0]*H + B[1]*F*H + C[2]*tx*F*H + C[3]*tx**2*F*H + C[4]*tx**3*F*H + C[5]*tx**4*F*H + C[6]*tx**5*F*H+\
               C[7]*tx*ycs*H + C[8]*tx**2*ycs*H + C[9]*tx**3*ycs*H + C[10]*tx**4*ycs*H + C[11]*tx**5*ycs*H + C[12]*ycs*H+\
               C[13]*F*ycs*H + C[14]*tx*F*ycs*H + C[15]*tx**2*F*ycs*H + C[16] + C[17]*tx*ycs + C[18]*tx**2*ycs +\
               C[19]*tx**3*ycs + C[20]*F*ycs + C[21]*tx*F*ycs + C[22]*tx**2*F*ycs
        DTC = sum_

    if ZHT > 240 and ZHT <= 300.0:
        H = 0.8
        sum_ = C[0]*H + B[1]*F*H + C[2]*tx*F*H + C[3]*tx**2*F*H + C[4]*tx**3*F*H + C[5]*tx**4*F*H + C[6]*tx**5*F*H +\
               C[7]*tx*ycs*H + C[8]*tx**2*ycs*H + C[9]*tx**3*ycs*H + C[10]*tx**4*ycs*H + C[11]*tx**5*ycs*H + C[12]*ycs*H+\
               C[13]*F*ycs*H + C[14]*tx*F*ycs*H + C[15]*tx**2*F*ycs*H + C[16] + C[17]*tx*ycs + C[18]*tx**2*ycs +\
               C[19]*tx**3*ycs + C[20]*F*ycs + C[21]*tx*F*ycs + C[22]*tx**2*F*ycs
        AA = sum_
        BB = C[0] + B[1]*F + C[2]*tx*F + C[3]*tx**2*F + C[4]*tx**3*F + C[5]*tx**4*F + C[6]*tx**5*F +\
             C[7]*tx*ycs + C[8]*tx**2*ycs + C[9]*tx**3*ycs + C[10]*tx**4*ycs + C[11]*tx**5*ycs + C[12]*ycs +\
             C[13]*F*ycs + C[14]*tx*F*ycs + C[15]*tx**2*F*ycs
        H = 3
        sum_ = B[0] + B[1]*F + B[2]*tx*F + B[3]*tx**2*F + B[4]*tx**3*F + B[5]*tx**4*F + B[6]*tx**5*F + \
               B[7]*tx*ycs + B[8]*tx**2*ycs + B[9]*tx**3*ycs + B[10]*tx**4*ycs + B[11]*tx**5*ycs + B[12]*H*ycs +\
               B[13]*tx*H*ycs + B[14]*tx**2*H*ycs + B[15]*tx**3*H*ycs + B[16]*tx**4*H*ycs + B[17]*tx**5*H*ycs + B[18]*ycs
        DTC300 = sum_
        sum_ = B[12]*ycs + B[13]*tx*ycs + B[14]*tx**2*ycs + B[15]*tx**3*ycs + B[16]*tx**4*ycs + B[17]*tx**5*ycs
        DTC300DZ = sum_
        CC = 3*DTC300 - DTC300DZ - 3*AA - 2*BB
        DD = DTC300 - AA - BB - CC
        ZP  = (ZHT-240)/60
        DTC = AA + BB*ZP + CC*ZP*ZP + DD*ZP*ZP*ZP

    if ZHT > 300 and ZHT <= 600:
        H = ZHT/100
        sum_ = B[0] + B[1]*F + B[2]*tx*F + B[3]*tx**2*F + B[4]*tx**3*F + B[5]*tx**4*F + B[6]*tx**5*F +\
               B[7]*tx*ycs + B[8]*tx**2*ycs + B[9]*tx**3*ycs + B[10]*tx**4*ycs + B[11]*tx**5*ycs + B[12]*H*ycs +\
               B[13]*tx*H*ycs + B[14]*tx**2*H*ycs + B[15]*tx**3*H*ycs + B[16]*tx**4*H*ycs + B[17]*tx**5*H*ycs + B[18]*ycs
        DTC = sum_

    if ZHT > 600 and ZHT <= 800.0:
        ZP = (ZHT - 600)/100
        HP = 6
        AA = B[0] + B[1]*F + B[2]*tx*F + B[3]*tx**2*F + B[4]*tx**3*F + B[5]*tx**4*F + B[6]*tx**5*F +\
             B[7]*tx*ycs + B[8]*tx**2*ycs + B[9]*tx**3*ycs + B[10]*tx**4*ycs + B[11]*tx**5*ycs + B[12]*HP*ycs +\
             B[13]*tx*HP*ycs + B[14]*tx**2*HP*ycs+ B[15]*tx**3*HP*ycs + B[16]*tx**4*HP*ycs + B[17]*tx**5*HP*ycs + B[18]*ycs
        BB  = B[12]*ycs + B[13]*tx*ycs + B[14]*tx**2*ycs + B[15]*tx**3*ycs + B[16]*tx**4*ycs + B[17]*tx**5*ycs
        CC  = -(3*AA+4*BB)/4
        DD  = (AA+BB)/4
        DTC = AA + BB*ZP + CC*ZP*ZP + DD*ZP*ZP*ZP
    return DTC

@jit(nopython=True)
def SEMIAN08 (DAY,HT,F10B,S10B,M10B):
    '''
    COMPUTE SEMIANNUAL VARIATION (DELTA LOG RHO)

    INPUT DAY, HEIGHT, F10B, S10B, M10B FSMB
          025.  650.   150.  148.  147. 151.
    OUTPUT FUNCTIONS FZ, GT, AND DEL LOG RHO VALUE

    DAY     (I)   DAY OF YEAR
    HT      (I)   HEIGHT (KM)
    F10B    (I)   AVE 81-DAY CENTERED F10
    S10B    (I)   AVE 81-DAY CENTERED S10
    M10B   (I)   AVE 81-DAY CENTERED M10
    FZZ     (O)   SEMIANNUAL AMPLITUDE
    GTZ     (O)   SEMIANNUAL PHASE FUNCTION
    DRLOG   (O)   DELTA LOG RHO
    '''
    TWOPI = Const.twopi

    # FZ GLOBAL MODEL VALUES
    # 1997-2006 FIT:
    FZM = np.array([0.2689,-0.01176, 0.02782,-0.02782, 0.3470e-3])

    # GT GLOBAL MODEL VALUES
    # 1997-2006 FIT:
    GTM = np.array([-0.3633, 0.08506, 0.2401,-0.1897, -0.2554,-0.01790, 0.5650e-3,-0.6407e-3,-0.3418e-2,-0.1252e-2])

    # COMPUTE NEW 81-DAY CENTERED SOLAR INDEX FOR FZ
    FSMB = F10B - 0.7*S10B - 0.04*M10B
    HTZ = HT/1e3
    FZZ = FZM[0] + FZM[1]*FSMB + FZM[2]*FSMB*HTZ + FZM[3]*FSMB*HTZ**2 + FZM[4]*FSMB**2*HTZ

    # COMPUTE DAILY 81-DAY CENTERED SOLAR INDEX FOR GT
    FSMB = F10B - 0.75*S10B - 0.37*M10B

    TAU = (DAY-1)/365
    SIN1P = np.sin(TWOPI*TAU)
    COS1P = np.cos(TWOPI*TAU)
    SIN2P = np.sin(2*TWOPI*TAU)
    COS2P = np.cos(2*TWOPI*TAU)

    GTZ = GTM[0] + GTM[1]*SIN1P + GTM[2]*COS1P + GTM[3]*SIN2P + GTM[4]*COS2P + GTM[5]*FSMB + GTM[6]*FSMB*SIN1P + GTM[7]*FSMB*COS1P + GTM[8]*FSMB*SIN2P + GTM[9]*FSMB*COS2P
    if FZZ < 1e-6: FZZ = 1e-6
    DRLOG = FZZ*GTZ

    return FZZ,GTZ,DRLOG