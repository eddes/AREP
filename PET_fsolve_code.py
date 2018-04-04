# PET calculation routine (2017)
# Edouard Walther (AREP, France) and Quentin Goestchel (ENS Paris-Saclay, France)
# based on: Peter Hoeppe's PET fortran code, from the VDI Norm 3787, Blatt 2
# and on : Djordje Spasic's python code downloaded on github:
# https://github.com/stgeorges/ladybug/commit/b0c2ea970252b62d22bf0e35d739db7f385a3f26
# the reference publication can be found on ResearchGate https://www.researchgate.net/publication/324168880_The_PET_comfort_index_Questioning_the_model
# (also available on Elsevier https://www.sciencedirect.com/science/article/pii/S0360132318301896 )


import numpy as np
import math as math
import scipy.optimize as optimize

# Skin and core temperatures set values
tc_set=36.6 # 36.8
tsk_set=34 # 33.7
tbody_set=0.1*tsk_set+0.9*tc_set # Calculation of the body temperature through a weighted average

# Skin blood flow calculation function:
def vasoC(tcore,tsk):
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
    # in the transient model, alpha is used to update tbody
    #alpha = 0.04177 + 0.74518 / (qmblood + 0.585417)
    alpha = 0.1
    return (qmblood,alpha)


# Sweating calculation function
def Suda(tbody,tsk):
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

# Vectorial MEMI balance calculation function
def Syst(T, Ta, Tmrt, HR, v, age, sex, ht, mbody, pos, M, icl,mode):
    # Conversion of T vector in an array
    arr = np.ones((3,1))
    arr[0,0]=T[0] #Corresponds to T_core
    arr[1,0]=T[1] #Corresponds to T_skin
    arr[2,0]=T[2] #Corresponds to T_clothes
    T=arr
    enbal_vec = np.zeros((3,1)) #required for the vectorial expression of the balance

    # Area parameters of the body:
    Adu = 0.203 * mbody ** 0.425 * ht ** 0.725
    feff=0.725
    if pos == 1 or pos == 3:
        feff = 0.725
    if pos == 2:
        feff = 0.696
    # Calculation of the Burton surface increase coefficient, k = 0.31 for Hoeppe:
    fcl = 1 + (0.31 * icl) # Increase heat exchange surface depending on clothing level
    facl = (173.51 * icl - 2.36 - 100.76 * icl * icl + 19.28 * icl ** 3.0) / 100
    Aclo = Adu * facl + Adu * (fcl - 1.0)
    Aeffr = Adu * feff  # Effective radiative area depending on the position of the subject
    
	# Partial pressure of water in the air depending on relative humidity and air temperature:
    if mode: # mode=True is the calculation of the actual environment
        vpa = HR / 100.0 * 6.105 * math.exp(17.27 * Ta / (237.7 + Ta )) #[hPa]
    else: # mode=False means we are calculating the PET
        vpa= 12 # [hPa] vapour pressure of the standard environment

    # Convection coefficient depending on wind velocity and subject position
    hc = 0
    if pos == 1:
        hc = 2.67 + (6.5 *v**0.67)
    if pos == 2:
        hc = 2.26 + (7.42 *v**0.67)
    if pos == 3:
        hc = 8.6 * (v ** 0.513)
    # modification of hc with the total pressure
	hc = hc * (p / po) ** 0.55

    # Base metabolism for men and women in [W]
    metab_female=3.19*mbody**0.75*(1.0+0.004*(30.0-age)+0.018*(ht*100.0/mbody**(1.0/3.0)- 42.1))
    metab_male=3.45*mbody**0.75*(1.0+0.004*(30.0-age)+0.01*(ht*100.0/mbody**(1.0/3.0)-43.4))
    # Source term : metabolic activity
    if mode==True: # = actual environment
		metab = (M + metab_male)/Adu
		fec = (M + metab_female)/Adu
    else:# False=reference environment
        metab = (80 + metab_male)/Adu
        fec = (80 + metab_female)/Adu
	
    he = 0.0
    # Attribution of internal energy depending on the sex of the subject
    if sex == 1:
        he = metab
    elif sex == 2:
        he = fec
    h = he *(1.0 - eta) # [W/m2]

    # Respiratory energy losses
    # Expired air temperature calculation:
    texp = 0.47 * Ta + 21.0  # [degC]
    # Pulmonary flow rate
    dventpulm = he * 1.44 * 10.0**(-6.0)
    # Sensible heat energy loss:
    eres = cair * (Ta - texp) * dventpulm  # [W/m2]
    # Latent heat energy loss:
    vpexp = 6.11 * 10.0**(7.45 * texp / (235.0 + texp))
    erel = 0.623 * Lvap / p * (vpa-vpexp) * dventpulm  # [W/m2]
    ere = eres + erel  # [W/m2]

    # Clothed fraction of the body approximation
    rcl = icl / 6.45  # Conversion in m2.K/W
    y = 0
    if facl > 1.0:
        facl = 1.0
    if icl >= 2.0:
        y = 1.0
    if icl > 0.6 and icl < 2.0:
        y = (ht - 0.2)/ht
    if icl <= 0.6 and icl > 0.3:
        y = 0.5
    if icl <= 0.3 and icl > 0.0:
        y = 0.1
    # calculation of the closing radius depending on the clothing level (6.28 = 2* pi !)
    r2 = Adu * (fcl - 1.0 + facl) / (6.28 * ht * y)  # External radius
    r1 = facl * Adu /(6.28 * ht * y)  # Internal radius
    di = r2 - r1
    # Calculation of the equivalent thermal resistance of body tissues
    alpha = vasoC(T[0,0],T[1,0])[1]
    tbody = alpha * T[1,0] + (1 - alpha) * T[0,0]
    htcl = (6.28 * ht * y * di)/(rcl * math.log(r2/r1)*Aclo)  # [W/(m2.K)]
    # Calculation of sweat losses
    qmsw = Suda(tbody,T[1,0])
    # Lvap/1000 = 2400 000[J/kg] divided by 1000 = [J/g] // qwsw/3600 for [g/m2/h] to [g/m2/s]
    esw = Lvap/1000* qmsw/3600  # [W/m2]
    # Saturation vapor pressure at temperature Tsk
    Pvsk = 6.105*math.exp((17.27 * (T[1,0]+273.15) - 4717.03)/ (237.7 + T[1,0])) # [hPa]
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
    i_m=0.38 # Woodcock's ratio
    R_ecl=(1/(fcl*hc) + rcl)/(Lw*i_m) # clothing vapour transfer resistance after Woodcock's method
    #R_ecl=0.79*1e7 # Hoeppe's method for E_diff
    ediff = (1 - w)*(Pvsk - vpa)/R_ecl  # diffusion heat transfer
    evap = -(ediff + esw)  # [W/m2]

    # Radiation losses
    # For bare skin area:
    rbare = Aeffr*(1.0 - facl) * emsk * sigm * ((Tmrt + 273.15)**(4.0) - (T[1,0] + 273.15)**(4.0))/Adu
    # For dressed area:
    rclo = feff * Aclo * emcl * sigm * ((Tmrt + 273.15)**(4.0) - (T[2,0] + 273.15)**(4.0))/Adu
    rsum = rclo+rbare

    # Convection losses #
    cbare = hc * (Ta - T[1,0]) * Adu * (1.0 - facl)/Adu  # [w/m^2]
    cclo = hc * (Ta - T[2,0]) * Aclo/Adu  # [W/m^2]
    csum = cclo+cbare

    # Balance equations of the 3-nodes model
    enbal_vec[0,0] = h + ere - (vasoC(T[0,0],T[1,0])[0]/3600*cb+5.28)*(T[0,0]-T[1,0]) # Core balance [W/m^2]
    enbal_vec[1,0] = rbare + cbare + evap + (vasoC(T[0,0],T[1,0])[0]/3600*cb+5.28)*(T[0,0]-T[1,0]) - htcl*(T[1,0]-T[2,0])  # Skin balance [W/m^2]
    enbal_vec[2,0] = cclo + rclo + htcl*(T[1,0]-T[2,0]) # Clothes balance [W/m^2]
    enbal_scal = h + ere + rsum + csum +evap

	#returning either the calculated core,skin,clo temperatures or the PET
    if mode:
		# if we solve for the system we need to return 3 temperatures 
        return [enbal_vec[0,0],enbal_vec[1,0],enbal_vec[2,0]]
    else:
		# solving for the PET requires the scalar balance only
        return enbal_scal

# Solving the 3 equation non-linear system
def resolution(Ta, Tmrt, HR, v, age, sex, ht, mbody, pos, M, icl, Tx):
    Tn = optimize.fsolve(Syst,Tx ,args=(Ta, Tmrt, HR, v, age, sex, ht, mbody, pos, M, icl, True))
    return (Tn, 1)

# PET calculation with dichotomy method 
def PET (age, sex, ht, mbody, pos, M, icl, Tstable,a,b,eps):
    # Definition of a function with the input variables of the PET reference situation
    def f(Tx):
        return Syst(Tstable, Tx, Tx, 50, 0.1, age, sex, ht, mbody, pos, M, 0.9,False)
    Ti = Tmin # Start of the search interval
    Tf = Tmax # End 	of the search interval
    pet = 0
    while Tf-Ti>eps: # Dichotomy loop
        if f(Ti)*f(pet)<0:
            Tf = pet
        else:
            Ti = pet
        pet = (Ti + Tf) / 2.0
    return pet

# Input data
# definition of constants
po = 1013.25 #[hPa]
rob = 1.06 # Blood density [kg/L]
cb = 3.64 * 1000. # Blood specific heat [J/kg/k]
cair = 1.01 * 1000. # Air specific heat  [J/kg/K]
emsk = 0.99 # Skin emissivity
emcl = 0.95 # Clothing emissivity
Lvap = 2.42 * 10. ** 6. # Latent heat of evaporation [J/Kg]
sigm = 5.67 * 10. ** (-8.) # Stefan-Boltzmann constant [W/(m2*K^(-4))]
eta = 0. # Body efficiency

# Initialisation of Temperature vector with respectively: Tcore, Tskin, Tcl
T = [38,40,40]
eps = 10**(-6) # numerical tolerance
# Dichotomy search interval (a=min / b=max)
Tmin = -40
Tmax = 60

# Input data for the PET 
Ta=30 # Air temperature in [oC]
Tmrt=60 # Mean radiant temperature in [oC]
HR=50 # Air relative humidity [%]
v=1 # Wind velocity [m/s]
age = 35
sex = 1 # 1 for men and 2 for women
pos = 1
mbody = 75 #[kg]
ht = 1.80 #[m]
p = 1013.25 #[hPa]
M = 80 # [W] Metabolic activity level
icl = 0.5 # [clo] Clothing level

# Results 
Tstable = resolution(Ta,Tmrt,HR,v,age,sex,ht,mbody,pos,M,icl,T)[0]
print("Nodes temperature [T_core, T_skin, T_clo]",Tstable)
print('Thermal Balance', Syst(Tstable, Ta, Tmrt, HR, v, age, sex, ht, mbody, pos, M, icl,True)[0])
print('PET:', round(PET(age, sex, ht, mbody, pos, M, icl, Tstable, Tmin, Tmax, eps),2))
