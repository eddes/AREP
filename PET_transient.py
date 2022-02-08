# -*- coding: utf-8 -*-
# Transient PET calculation routine (work in progress as of Feb. 2022)
# by Edouard Walther and William Dumontaud 2021-2022
# based on the paper "The PET comfort index: Questioning the model"
# E. Walther and Q. Goetschel 
# in Building & Environment (2018)

import numpy as np
import math as math
from scipy.optimize import fsolve

# fonction for vapour pressure
def fc_pvs(T):
    """
    Vapour pressure

    Parameters
    ----------
    T : float
        Dry bulb temperature.

    Returns
    -------
    pvs : float
        vapour pressure at saturation.
    """
    a,b,c,d=0.07252,0.0002881,0.00000079,611
    pvs=d*np.exp(a*T -b*np.power(T,2)+ c*np.power(T,3))
    return pvs


########### Useful functions for code ##########

# - Function A_dubois
def A_dubois(mbody,ht):
    """
    A function to compute the Dubois body surface area.

    Parameters
    ----------
    mbody : float
        body mass in kg.
    ht : float
        body height in meters.

    Returns
    -------
    float 
        the Dubois (naked) body area.

    """
    return 0.203 * mbody ** 0.425 * ht ** 0.725

# - Skin blood flow calculation function:

def vasoC(tcore,tsk):
    """
    Defines the vasomotricity (blood flow) in function of the core and skin temperatures.
    
    Parameters
    ----------
    tcore : float
        The body core temperature.
    tsk : float
        The body skin temperature.

    Returns
    -------
    qmblood : float
        Blood flow rate (kg/m2/h).
    alpha : float
        Repartition of body mass between core and skin (-).

    """
    # skin and core temperatures set values
    tc_set=36.6 # 36.8
    tsk_set=34 # 33.7
    #Set value signals
    sig_skin = tsk_set - tsk
    sig_core = tcore - tc_set
    if sig_core<0:
        # In this case, Tcore<Tc_set --> the blood flow is reduced
        sig_core=0.
    if sig_skin<0:
        # In this case, Tsk>Tsk_set --> the blood flow is increased
        sig_skin=0.
    # 6.3 L/m^2/h is the set value of the blood flow
    qmblood = (6.3 + 75. * sig_core) / (1. + 0.5 * sig_skin)
    # 90 L/m^2/h is the blood flow upper limit
    if qmblood>90:
        qmblood=90.
    # in other models, alpha is used to update tbody
    alpha = (0.04177 + 0.74518) / (qmblood + 0.585417)
    #alpha = 0.1
    return (qmblood,alpha)

# - Sweating calculation function:
def Suda(tbody,tsk):
    """
    Defines the sweating mechanism depending on the body and core temperatures.
    
    Parameters
    ----------
    tbody : float
        The actual body temperature (weighted average between skin and core temperatures).
    tsk : float
        The actual skin temperature.

    Returns
    -------
    qmsw : float
        The sweating flow rate in g/m2/h.

    """
    tc_set=36.6 # 36.8
    tsk_set=34 # 33.7
    tbody_set=0.1*tsk_set+0.9*tc_set # Calculation of the body temperature through a weighted average
    sig_body = tbody - tbody_set
    sig_skin = tsk - tsk_set
    if sig_body<0:
        #In this case, Tbody<Tbody_set --> The sweat flow is 0
        sig_body=0.
    if sig_skin<0:
        # In this case, Tsk<Tsk_set --> the sweat flow is reduced
        sig_skin=0.
    #qmsw = 170 * sig_body * math.exp((sig_skin) / 10.7)  # [g/m2/h] is the expression from Gagge's model
    qmsw = 304.94*10**(-3) * sig_body
    # 500 g/m^2/h is the upper sweat rate limit
    if qmsw > 500:
        qmsw = 500
    return (qmsw)

# - Function for convective transfert coefficient depending
def fc_hc(pos,v,p):
    """
    Computes the convective heat transfer coefficient around the body.
    
    Parameters
    ----------
    pos : int
        position of the inividual.
    v : float
        Air velocity magnitude.
    p : TYPE
        total atmospheric pressure.

    Returns
    -------
    float
        hc the convective heat transfer coefficient.
    """
    hc = 0
    if pos==1: # sitting
        hc = 2.67 + 6.5*v**0.67
    if pos==2: # standing
        hc = 2.26 + 7.42*v**0.67
    if pos==3: # crouching
        hc = 8.6*v**0.513
    # modification of hc with the total pressure
    return hc*(p/po)**0.55

# - Function calculating the clothed area
def fc_Aclo(icl,Adu):
    """
    Defines the clothed body surface area depending on the clothing insulation level.
    
    Parameters
    ----------
    icl : TYPE
        The clothing insulation in clo .
    Adu : TYPE
        Dubois body surface area.

    Returns
    -------
    Aclo : TYPE
        The clothed surface area for external heat transfer.

    """
    fcl = 1 + (0.31 * icl) # Increase heat exchange surface depending on clothing level
    facl = (173.51 * icl - 2.36 - 100.76 * icl * icl + 19.28 * icl ** 3.0) / 100
    Aclo = Adu * facl + Adu * (fcl - 1.0)
    return Aclo

# - Function allowing to obtain the effective radiative coefficient:
def fc_feff(pos):
    """
    Returns the effective radiation coefficient of individuals
    
    Parameters
    ----------
    pos : int
        Position of the individual (1=sitting, 2=standing, 3=crouching)

    Returns
    -------
    feff : TYPE
        The effective radiation coefficient of individuals.
        
    """
    feff=0.725
    if pos == 1 or pos == 3:
        feff = 0.725
    if pos == 2:
        feff = 0.696
    return feff

# - Function to obtain the equivalent conductance coefficient between the skin and the clothing:
def fc_htcl(icl,ht,Adu,Aclo):
    """
    Gives the equivalent conductance between skin and clothing.
    
    Parameters
    ----------
    icl : float
        the clothing level of insulation.
    ht : float
        individual height .
    Adu : float 
        The Dubois surface area.
    Aclo : float
        The clothing surface area.

    Returns
    -------
    htcl : float
        The equivalent conductance between skin and clothing.

    """
    fcl = 1 + (0.31 * icl) # Increase heat exchange surface depending on clothing level
    facl = (173.51 * icl - 2.36 - 100.76 * icl * icl + 19.28 * icl ** 3.0) / 100
    # Clothed fraction of the body approximation
    rcl = icl / 6.45  # Conversion in m2.K/W
    y = 0
    if facl > 1.0: facl = 1.0
    if icl >= 2.0: y = 1.0
    if icl > 0.6 and icl < 2.0: y = (ht - 0.2)/ht
    if icl <= 0.6 and icl > 0.3: y = 0.5
    if icl <= 0.3 and icl > 0.0: y = 0.1
    # calculation of the clothing radius depending on the clothing level (6.28 = 2* pi !)
    r2 = Adu * (fcl - 1.0 + facl) / (6.28 * ht * y)  # External radius
    r1 = facl * Adu /(6.28 * ht * y)  # Internal radius
    di = r2 - r1
    # Calculation of the equivalent thermal resistance of body tissues
    htcl = (6.28 * ht * y * di)/(rcl * math.log(r2/r1)*Aclo)  # [W/(m2.K)]
    return htcl

# - The clothing model of Schiavon and Lee (clo=clo(T 6h du mat))
def calcul_clo_T6h(T_6AM):
    """
    A function to obtain the clo level from the 6AM outdoor air temperature after Schiavon & Lee's correlation.
    
    Parameters
    ----------
    T_6AM : float
        Temperature at 6AM on the morning of the day the comfort computation if performed.

    Returns
    -------
    clo : float
        clo level after Schiavon & Lee's correlation.

    """
    if T_6AM > 5. and T_6AM < 26.:
        clo = 10**(-0.1635 - 0.0066 * T_6AM)
    elif T_6AM > -5 and T_6AM <= 5:
        clo = 0.818 - 0.0364 * T_6AM
    elif T_6AM >= 26:
        clo = 0.46
    elif T_6AM <= -5:
        clo = 1.
    return clo


def fc_clo_properties(v_walk,v,he,icl,i_m):
    """
    An integration of the correction of the clothin heat and mass transfer properties after Havenith & Holmer 2002.

    Parameters
    ----------
    v_walk : float
        Walking velocity.
    v : float
        Air velocity.
    he : float
        Total metabolic activity (W/m2).
    icl : float
        Clothing insulation level.
    i_m : float
        Woodcock ratio for mass transfer.

    Returns
    -------
    v : float
        Effective air velocity.
    icl : float
        Corrected clothing insulation level.
    i_m : float
        Corrected Woodcock ratio.
    """
    #correction of clo properties after (Havenith et al. 2002) & (Holmer et.al 2002)
    # v_walk depending on metabolic rate
    if v_walk < 0.7:
        v_walk = 0.0052 * (he/58 - 58)
    # "effective" wind speed including body movement
    v = math.sqrt(v**2 + v_walk**2)
    # correction 
    corr_T = math.exp(0.042 - 0.398 * v  + 0.066 * v**2 - 0.378 * v_walk + 0.094 * v_walk**2)
    # upper bound for correlation validity
    if v > 3.5 or corr_T>0.582:
        corr_T = 0.582
    icl = icl * corr_T  # corrected insulation level
    i_m = i_m* (4.9 - 6.5 * corr_T + 2.6 * corr_T**2) # corrected Woodcock ratio
    return v,icl,i_m 

########## Vectorial MEMI balance calculation function ##########
def Syst(T, Ta, Tmrt, HR, v, age, sex, ht, mbody, pos, M, icl, mode, p):
    """
    This function allows to solve for the PET : either it solves the vectorial balance of the 3 unknown temperatures (Tcore, Tskin, Tclothing) or it solves for the environment operative temperature that would yield the same energy balance as the actual environment.
    
    Parameters
    ----------
    T : float
        [Tcore, Tskin, Tclothing].
    Ta : float
        Dry bulb air temperature.
    Tmrt : float
        Mean radiant temperature.
    HR : float
        Relative humidity (0-100%).
    v : float
        Air velocity (m/s).
    age : float
        age.
    sex : int
        male (1) or female (2).
    ht : float
        body height (m).
    mbody : float
        body mass (kg).
    pos : int
        position.
    M : float
        Metabolic heat (W).
    icl : float
        clothing level (clo).
    mode : boolean
        True=solve 3eqs/3unknowns, False=solve for PET.
    p : float
        Atmospheric pressure (hPa).

    Returns
    -------
    float
        PET or energy balance.
    """
    T=np.reshape(T,(3,1)) # reshape to proper dimensions for fsolve
    enbal_vec = np.zeros((3,1)) #required for the vectorial expression of the balance
    # Area parameters of the body:
    Adu = A_dubois(mbody,ht)
    # Base metabolism for men and women in [W]
    metab_female=3.19*mbody**0.75*(1.0+0.004*(30.0-age)+0.018*(ht*100.0/mbody**(1.0/3.0)- 42.1))
    metab_male=3.45*mbody**0.75*(1.0+0.004*(30.0-age)+0.01*(ht*100.0/mbody**(1.0/3.0)-43.4))
    # Source term : metabolic activity
    if mode==True: # = actual environment
        metab_M = (M + metab_male)/Adu
        metab_F = (M + metab_female)/Adu
    else:# False=reference environment
        metab_M = (80 + metab_male)/Adu
        metab_F = (80 + metab_female)/Adu
    he = 0.0 # total internal energy
    # Attribution of internal energy depending on the sex of the subject
    if sex == 1: he = metab_M
    else: he = metab_F
    # impact of efficiency
    h = he *(1.0 - eta) # [W/m2]

    # correction for wind 
    i_m = 0.38 # Woodcock ratio for vapour transfer through clothing [-]
    v_walk=0 # walking velocity
    v, icl, i_m = fc_clo_properties(v_walk, v, he, icl, i_m)
    
    # Calculation of the Burton surface increase coefficient, k = 0.31 for Hoeppe:
    fcl = 1 + (0.31 * icl) # Increase heat exchange surface depending on clothing level
    facl = (173.51 * icl - 2.36 - 100.76 * icl * icl + 19.28 * icl ** 3.0) / 100
    Aclo = fc_Aclo(icl,Adu)
    feff = fc_feff(pos) # effective radiation factor
    Aeffr = Adu*feff  # Effective radiative area depending on the position of the subject

    # Partial pressure of water in the air depending on relative humidity and air temperature:
    if mode: # mode=True is the calculation of the actual environment
        vpa = HR/100.0 * fc_pvs(Ta)/100 #[hPa]
    else: # mode=False means we are calculating the PET
        vpa= 12 # [hPa] vapour pressure of the standard environment

    # Convection coefficient depending on wind velocity and subject position
    hc = fc_hc(pos,v,p)

    # Respiratory energy losses
    # Expired air temperature calculation:
    texp = 0.47 * Ta + 21.0  # [degC]
    # breathing flow rate
    dventpulm = he * 1.44 * 10.0**(-6.0)
    # Sensible heat energy loss:
    eres = cair * (Ta - texp) * dventpulm  # [W/m2]
    # Latent heat energy loss:
    vpexp=fc_pvs(texp)/100 # hPa
    erel = 0.623 * Lvap / p * (vpa-vpexp) * dventpulm  # [W/m2]
    ere = eres + erel  # [W/m2]
    
    # Clothed fraction of the body approximation
    rcl = icl / 6.45  # Conversion in m2.K/W
    # Calculation of the equivalent thermal resistance of body tissues
    alpha = vasoC(T[0,0],T[1,0])[1]
    tbody = alpha * T[1,0] + (1 - alpha) * T[0,0]
    htcl = fc_htcl(icl,ht,Adu,Aclo)  # [W/(m2.K)]
    # Calculation of sweat losses
    qmsw = Suda(tbody,T[1,0])
    # Lvap/1000 = 2400 000[J/kg] divided by 1000 = [J/g] // qwsw/3600 for [g/m2/h] to [g/m2/s]
    esw = Lvap/1000* qmsw/3600  # [W/m2]
    # Saturation vapor pressure at temperature Tsk
    Pvsk =fc_pvs(T[1,0])/100 # hPa
    # Calculation of vapour transfer
    Lw = 16.7*10 **(-1)  # [K/hPa] Lewis factor
    he_diff = hc * Lw # diffusion coefficient of air layer
    fecl=1/(1+0.92*hc*rcl) # Burton efficiency factor
    emax = he_diff * fecl * (Pvsk - vpa) # maximum diffusion at skin surface
    w = esw / emax  # skin wettedness
    if w > 1:
        w=1
        delta = esw-emax
        if delta < 0:
            esw=emax
    if esw < 0:
        esw=0
    #i_m= Woodcock's ratio (see above)
    R_ecl=(1/(fcl*hc) + rcl)/(Lw*i_m) # clothing vapour transfer resistance after Woodcock's method
    #R_ecl=0.79*1e7 # Hoeppe's method for E_diff (unreferenced and rather weird...)
    ediff = (1 - w)*(Pvsk - vpa)/R_ecl  # diffusion heat transfer
    evap = -(ediff + esw)  # [W/m2]

    # Radiation losses
    #... for bare skin
    rbare = Aeffr*(1.0 - facl) * emsk * sigm * ((Tmrt + 273.15)**(4.0) - (T[1,0] + 273.15)**(4.0))/Adu
    #... for clothed area
    rclo = feff * Aclo * emcl * sigm * ((Tmrt + 273.15)**(4.0) - (T[2,0] + 273.15)**(4.0))/Adu
    rsum = rclo+rbare # radiation total

    # Convection losses
    #... for bare skin
    cbare = hc*(Ta - T[1,0]) * Adu * (1.0 - facl)/Adu  # [W/m^2]
    #... for clothed area
    cclo =  hc*(Ta - T[2,0]) * Aclo/Adu  # [W/m^2]
    csum = cclo+cbare # convection total

    # Balance equations of the 3-nodes model
    enbal_vec[0,0] = h + ere - (vasoC(T[0,0],T[1,0])[0]/3600*cb+5.28)*(T[0,0]-T[1,0]) # Core balance [W/m^2]
    enbal_vec[1,0] = rbare + cbare + evap + (vasoC(T[0,0],T[1,0])[0]/3600*cb+5.28)*(T[0,0]-T[1,0]) - htcl*(T[1,0]-T[2,0])  # Skin balance [W/m^2]
    enbal_vec[2,0] = cclo + rclo + htcl*(T[1,0]-T[2,0]) # Clothes balance [W/m^2]
    enbal_scal = h + ere + rsum + csum + evap

    #returning either the calculated core,skin,clo temperatures or the PET
    if mode==True:
        # if we solve for the system we need to return 3 temperatures
        return [enbal_vec[0,0],enbal_vec[1,0],enbal_vec[2,0]]
    else:
        # solving for the PET requires the scalar balance only
        return enbal_scal

########## PET calculation with dichotomy method ##########
def fc_PET(age, sex, ht, mbody, pos, M, icl, Tstable,p):
    """
    This function solves for the PET using scipy's fsolve
    
    Parameters
    ----------
    age : float
        age of individual.
    sex : int
        sex of individual (1=male, 2=female).
    ht : float
        individual height (m).
    mbody : float
        individual mass (kg).
    pos : int
        body position (1,2,3=sitting, standing, crouching).
    M : float
        metabolic heat (W).
    icl : float
        clothing insulation level (clo).
    Tstable : float
        3 temperatures obtained from the actual environment (Tcore,Tskin,Tclo).
    p : float
        atmospheric pressure (hPa).

    Returns
    -------
    float
        The PET comfort index.

    """
    # Definition of a function with the input variables of the PET reference situation
    def f(Tx):
        return Syst(Tstable, Tx, Tx, 50, 0.1, age, sex, ht, mbody, pos, M, 0.9,False,p)
    # solving for PET 
    PET_guess=Tstable[2] # start with the clothing temperature
    PET=fsolve(f,PET_guess)
    return PET

########## Function calculating the PET in transient regime ##########
def fc_transient_PET(Tinit,Ta,Tmrt,HR,v,icl,p,Time,dt):
    """
    This function computes the PET in transient regime.

    Parameters
    ----------
    Tinit : float
        The initial set of temperatures: Tcore,Tskin,Tclo
    Ta : float
        Dry bulb air temperature.
    Tmrt : float
        Mean radiant temperature.
    HR : float
        Relative humidity (0-100%).
    v : float
        Air velocity (m/s).
    icl : float
        Clothing level (clo).
    p : float
        Atmospheric pressure.
    Time : float
        Exposure time of the individual to the ambient conditions (s).
    dt : float
        Simulation time step for the explicite Euler scheme (s). Expect oscillations for large dt and elevated temperatures (high sweating rates+evaporation causes instability).

    Returns
    -------
    TYPE
        PET at the end of the simulation and simulation time.

    """
    # Initialisaing the 3-nodes temperature array over time
    Temp = np.zeros((3,int(Time/dt+1)))
    # Initial values for Tcore,Tskin,Tclothing
    Temp[:,0] = Tinit
    
    # time for the simulation
    time=np.linspace(0,Time,int(Time/dt+1))
    # array that will contain the calculated PETs
    Temp_pet = np.zeros(len(time))
    
    # Initialisation for the PET
    Temp_pet[0] = fc_PET(age, sex, ht, mbody, pos, M, icl, Temp[:,0],p)
    
    S = A_dubois(mbody,ht)
    Aclo = fc_Aclo(icl,S)
    
    # a function to be able to 'fsolve' Tclo
    def fc_Tcl(Tcl,Tsk):
        a = fc_hc(pos,v,p)*Aclo/S
        b = Aclo*fc_feff(pos)*emcl*sigm/S
        htcl = fc_htcl(icl,ht,S,Aclo)
        return -Tcl + a/(a+htcl)*Ta + b/(a+htcl)*((Tmrt+273.15)**4-(Tcl+273.15)**4) + htcl/(a+htcl)*Tsk
    
    # time loop 
    for k in range(len(time)-1):
        # Evaluation de la température du noyau à l'instant t+1
        alphA = vasoC(Temp[0,k],Temp[1,k])[1]
        # solve for Tcore and Tskin
        TcoreTskin=Syst(Temp[:,k], Ta, Tmrt, HR, v, age, sex, ht, mbody, pos, M, icl, True ,p)
        # explicit Euler forward for the computation of Tcore
        Temp[0,k+1] = Temp[0,k] + dt*S/(mbody*cp*(1-alphA))*TcoreTskin[0]
        # explicit Euler forward for the computation of Tskin
        Temp[1,k+1] = Temp[1,k] + dt*S/(mbody*cp*alphA)*TcoreTskin[1]
        # Steady state for Tclothing w/ previous temperature as initial guess
        Temp[2,k+1] = fsolve(fc_Tcl, Temp[2,k], args=(Temp[1,k+1]))
        # out of the loop if only the resulting PET matters
        Temp_pet[k+1] = fc_PET(age, sex, ht, mbody, pos, M, icl, Temp[:,k+1], p)

    # return Temp_pet, time
    return Temp_pet[-1], time[-1]

########## Input data ##########

# Definition of constants:
po = 1013.25 #[hPa]
rob = 1.06 # Blood density [kg/L]
cb = 3640 # Blood specific heat [J/kg/k]
cair = 1010 # Air specific heat  [J/kg/K]
emsk = 0.99 # Skin emissivity
emcl = 0.95 # Clothing emissivity
Lvap = 2.42*10**6 # Latent heat of evaporation [J/Kg]
sigm = 5.67*10**(-8) # Stefan-Boltzmann constant [W/(m2*K^(-4))]
eta = 0. # Body efficiency

# Input data for the PET
age = 35
sex = 1 # 1 for men and 2 for women
pos = 1 # 1=standing
mbody = 75 # body mass [kg]
ht = 1.80 # body height [m]
M = 100 # [W] Metabolic activity level
cp=3492 # body specific heat capacity

envmnt="office"
# envmnt="high_MRT"
# envmnt="winter"

if envmnt=="office":
    Tair, Tmrt = 20, 20
    HR, v, p, T6h = 50, .15, 1013.25, 15
if envmnt=="high_MRT":
    Tair, Tmrt = 30, 60
    HR, v, p, T6h = 20, .5, 1018.25, 15
if envmnt=="winter":
    Tair, Tmrt = 0,0
    HR, v, p, T6h = 20, .5, 1018.25, 15

# get the clothing insulation according to the weather
icl = calcul_clo_T6h(T6h)
dt=10
Time=3600
Tinit = [36.8,33.8,(33.8+(Tair+Tmrt)/2)*0.5]

# compute 
Temp_PET, time = fc_transient_PET(Tinit,Tair,Tmrt,HR,v,icl,p,Time,dt)

# decide whether print or plot
if sum(np.shape(Temp_PET))==0:
    print(round(Temp_PET,2))
else:
    import matplotlib.pyplot as plt
    plt.plot(Temp_PET)
