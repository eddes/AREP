# -*- coding: utf-8 -*-
#quelques couleurs
rouge_A='#C60C2E'
vert1_A='#005157'
vert2_A='#627D77'
vert3_A='#9EB28F'
vert4_A='#C5E5A4'
gris1_A='#595A5C'
coule=[rouge_A,vert1_A,vert2_A,vert3_A,vert4_A,gris1_A]
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 21:12:00 2022

@author: walthere
"""
import numpy as np
from scipy.optimize import fsolve

def fc_clo_properties(v_air,v_walk,clo):
    """
    This function provides the clothing insulation levels after ISO 9920,
    i.e. the heat and evaporative transfer resistance of clothing depending on air velocity and walking.
    Range of validity : 
            0.15 < v_air < 3.5 m/s
            0 < v_walk < 1.2 m/s
            0.6 < clo < 2 clo
    
    Parameters
    ----------
    v_air : float
        air velocity magnitude (m/s).
    v_walk : float
        walking velocity (m/s).
    clo : float
        clothing insulation level in (clo).

    Returns
    -------
    Rclo_mod : float
        modified clothing insulation resistance (m2.K/W).
    Recl_mod : float
        modified clothing vapour resistance (m2.Pa/W).

    """
    # validity ranges for air and walking velocities
    v_air=max(0.15,v_air)
    v_walk=min(v_walk, 1.2)
    # incidence on HEAT TRANSFER
    dv=v_air-0.15 # m/s
    # ratio for clo # after ISO 9920
    k=np.exp(-0.281*dv+0.044*dv**2 -0.492*v_walk +0.176*v_walk**2)
    # ratio superficial heat transfer coefficient # after ISO 9920
    ka=np.exp( -0.533 *dv + 0.069*dv**2 + 0.201*v_walk**2)
    
    # 1 clo = 0.155 m2.K/W
    Rclo=clo*0.155 # m2.K/W
    Rclo_mod=k*Rclo
    # incidence on EVAPORATIVE MASS TRANSFER
    i_cl=0.34 # vapour permation index of clothing
    Lewis=0.0165 # lewis ratio K/Pa
    Recl = Rclo/(i_cl*Lewis) # Pa/K
    km = 0.3 -0.5*k +1.2*k**2
    Recl_mod = km*Recl
    
    sigma=5.67*1e-8 # Boltzmann's constant
    hc=fc_hc(v_air,1)
    Ra=1/(hc+sigma*0.73*0.95*293**3) # total conv/rad heat trasnfer resistance
    Ra_mod=ka*Ra # total modified surface resistance
    
    return Rclo_mod, Recl_mod

def fc_hc_properties(v_air,v_walk, hr):
    """
    This function provides the air lever insulation levels after ISO 9920,
    Range of validity : 
            0.15 < v_air < 3.5 m/s
            0 < v_walk < 1.2 m/s
            0.6 < clo < 2 clo
    
    Parameters
    ----------
    v_air : float
        air velocity magnitude (m/s).
    v_walk : float
        walking velocity (m/s).
    
    Returns
    -------
    Ra_mod : float
        modified superficial heat tranfer resistance (m2.K/W).

    """
    # validity ranges for air and walking velocities
    v_air=max(0.15,v_air)
    v_walk=min(v_walk, 1.2)
    # incidence on HEAT TRANSFER
    dv=v_air-0.15 # m/s
    # ratio superficial heat transfer coefficient # after ISO 9920
    ka=np.exp( -0.533 *dv + 0.069*dv**2 + 0.201*v_walk**2)
    
    hc=fc_hc(v_air,1)
    Ra=1/(hc+hr) # total conv/rad heat trasnfer resistance
    Ra_mod=ka*Ra # total modified surface resistance
    hc_mod=1/Ra_mod-hr
    return Ra_mod, hc_mod

# - The clothing model of Schiavon and Lee (clo=clo(T 6h du mat))
def fc_clo_T6am(T_6AM):
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


def w(pv, p):
    return 0.622 * pv/(p - pv)

def Cp_ah(pv, p):
    cpa = 1006
    cpv = 1830
    water = w(pv, p)
    return (cpa + water * cpv)/(1 + water)

def v_spe(T, pv, p):
    return (461.24 *(T+273.15)*(0.622 + w(pv, p)))/p

def fc_Lewis(T,pv,p):
    vs=v_spe(T, pv, p)
    Cp=Cp_ah(pv, p)
    Lewis = 2434*vs/(Cp*1.04*np.power(0.83,(2/3)))*(18/8.32/(T+ 273.15))
    return Lewis

def fc_pvs(T):
    pvs = np.exp(18.6686-4030.183/(T+ 235))*133.3223684  # conversion mmHg>Pa
    return pvs

def fc_hc(v,p_in_atm):
    # h_cc corrected convective heat transfer coefficient
    h_cc = 3.0*np.power(p_in_atm, 0.53)
    # h_fc forced convective heat transfer coefficient, W/(m2.K)
    h_fc = 8.600001*np.power((v* p_in_atm), 0.53)
    h_cc = max(h_cc, h_fc)
    return h_cc

def set_optimized(tdb,tr,HR,v,T6h,p_atmospheric, iso_9920):
    met=1
    v_walk=0
    # validation=True
    validation=False
       
    # numerical setup (time step, time and degree of "implicitness")
    dt,tmax,theta=120,3600,1
    # conversion en Pa
    vapor_pressure=HR/100*fc_pvs(tdb)
    if validation==False:
        # clo after Schiavon and Lee 2014
        i_m=0.45
        i_cl=0.38
        clo=fc_clo_T6am(T6h)
        # variation clo=clo(vair)
    else:
        i_cl=0.45
        i_m=0.45
        clo=0.5
        
    wme=0 # work mean efficiency
    body_surface_area=1.78
    body_position=1
    calculate_ce=False,
    max_skin_blood_flow=90
    
    # variables to check if person is experiencing heat strain
    heat_strain_blood_flow = False  # reached max blood flow
    heat_strain_sweating = False  # reached max regulatory sweating
    heat_strain_w = False  # reached max skin wettedness
    
    # cp pour le calcul de stockage
    cp_bl=4200 # J/kg/K
    cp_body=3492 # J/kg/K
    Lv=2426 # J/g_H2O latent heat of vapourisation
    cp_air=1006 # J/kg/K
    cp_vap=1830# J/kg/K
    # Initial variables as defined in the ASHRAE 55-2017
    air_speed = max(v, 0.1)
    k_clo = 0.25
    body_weight = 66  # body weight in kg
    met_factor = 58.2  # met conversion factor
    sbc = 5.67*1e-8  # Stefan-Boltzmann constant (W/m2K4)
    c_sw = 170  # driving coefficient for regulatory sweating
    c_dil = 200  # driving coefficient for vasodilation ashrae says 50 see page 195
    c_str = 0.5  # driving coefficient for vasoconstriction

    temp_skin_neutral = 33.7
    temp_core_neutral = 36.8
    alfa = 0.1
    temp_body_neutral = alfa*temp_skin_neutral + (1-alfa)*temp_core_neutral
    skin_blood_flow_neutral = 6.3 # L/m2/h
    
    t_skin = temp_skin_neutral
    t_core = temp_core_neutral
    m_bl = skin_blood_flow_neutral

    # initialize some variables
    e_skin = 0.1 * met  # total evaporative heat loss, W
    q_sensible = 0  # total sensible heat loss, W
    w = 0  # skin wettedness
    _set = 0  # standard effective temperature
    e_rsw = 0  # heat lost by vaporization sweat
    e_diff = 0  # vapor diffusion through skin
    e_max = 0  # maximum evaporative capacity
    m_rsw = 0  # regulatory sweating
    q_res = 0  # heat loss due to respiration

    pressure_in_atmospheres = p_atmospheric/101325
    time=np.arange(0,tmax,dt)

    r_clo = 0.155 * clo  # thermal resistance of clothing, C M^2 /W
    f_a_cl = 1.0 + 0.15* clo  # increase in body surface area due to clothing
    # lewis ratio scaled with pressure
    # generally about lr = 16.7*1e-3
    lr=fc_Lewis(tdb,vapor_pressure,p_atmospheric)/pressure_in_atmospheres
    rm = met * met_factor  # metabolic rate
    m =  met * met_factor  # metabolic rate

    if clo <= 0:
        w_max = 0.38*np.power(air_speed, -0.29)  # critical skin wettedness
        i_cl = 1.0  # permeation efficiency of water vapour through the clothing layer
    else:
        w_max = 0.59*np.power(air_speed, -0.08)  # critical skin wettedness
    
    h_cc = fc_hc(v, pressure_in_atmospheres)
    
    if not calculate_ce and met > 0.85:
        h_c_met = 5.66 * (met - 0.85)**0.39
        h_cc = max(h_cc, h_c_met)

    h_r = 4.7  # linearized radiative heat transfer coefficient
    h_t = h_r + h_cc  # sum of convective and radiant heat transfer coefficient W/(m2*K)
    r_a = 1.0 / (f_a_cl * h_t)  # resistance of air layer to dry heat
    t_op = (h_r * tr + h_cc * tdb) / h_t  # operative temperature
    h_c0 = 3.0*np.power(pressure_in_atmospheres, 0.53)
    
    # compute r_clo and r_eclo with the ISO 9920 method
    # mettre en argument
    #   => lewis
    #   => i_cl car \ne 0.38
    lr=fc_Lewis(tdb,vapor_pressure,p_atmospheric)/pressure_in_atmospheres 
    if iso_9920:
        # methode i_cl
        r_clo, r_eclo = fc_clo_properties(v, v_walk, clo)
        h_r = 4.85
        r_a, h_cc = fc_hc_properties(v, v_walk, h_r)
        r_ea = 1 /(lr * f_a_cl * h_cc)  # evaporative resistance air layer
        r_a = r_a/f_a_cl
    else:
        # methode i_cl
        r_ea = 1.0 / (lr * f_a_cl * h_cc)  # evaporative resistance air layer
        r_eclo = r_clo / (lr * i_cl)
        
    r_etot = (r_ea+ r_eclo)
    
    def fc_TcoreTskin(t_core, t_skin, qsens_return,k):
        #update hc
        air_speed = max(v, 0.1)
        r_a, h_cc =fc_hc_properties(air_speed , v_walk, h_r)
        # t_cl temperature of the outer surface of clothing
        t_cl = (r_a * t_skin + r_clo * t_op) / (r_a + r_clo)  # initial guess
        # built-in function for Tclo (consider solving the proper t_cl**4 equation instead)
        def fc_Tcl(t_cl):
            if body_position == "sitting": ratio=0.7
            else:ratio=0.73
            h_r = 4.0 * 0.95 * sbc * ((t_cl + tr)/2.0 + 273.15) ** 3.0 *ratio
            r_a, h_ccm = fc_hc_properties(air_speed, v_walk, h_r)
            h_t = h_r + h_ccm
            r_a = 1.0 / (f_a_cl * h_t)
            t_op = (h_r * tr + h_cc * tdb) / h_t            
            
            # h_t = h_r + h_cc
            # r_a = 1.0 / (f_a_cl * h_t)
            # t_op = (h_r * tr + h_cc * tdb) / h_t
            return - t_cl + (r_a * t_skin + r_clo * t_op) / (r_a + r_clo)
        r_a=r_a/f_a_cl
        #solve for Tclo
        t_cl=fsolve(fc_Tcl,t_cl)
        # fluxes
        q_sensible = (t_skin - t_op) / (r_a + r_clo)  # total sensible heat loss, W
        # hf_cs rate of energy transport between core and skin, W
        hf_cs = (t_core - t_skin) * (5.28 + cp_bl*m_bl/3600)
        # pour Tresp Thellier 1989
        t_resp=32.5 + 0.066*tdb + 1.99*1e-6*vapor_pressure 
        # respiration mass flow rate (2.14*1e-5)
        qm_resp=m*0.3/3600/58.2 # kg/m2/s
        # sensible respiration losses (stringent)
        wair=0.622*vapor_pressure/(p_atmospheric-vapor_pressure) # kgw/kga
        c_res = qm_resp*(cp_air+cp_vap*wair)*(t_resp-tdb)  # convective heat loss respiration
        # latent respiration losses (stringent)
        pv_resp=fc_pvs(t_resp)
        wresp=0.622*pv_resp/(p_atmospheric-pv_resp) # kgw/kga
        q_res = 0.622*qm_resp*Lv*(wresp - wair)  # heat loss due to respiration

        s_core = m - hf_cs - q_res - c_res - wme  # rate of energy storage in the core
        s_skin = hf_cs - q_sensible - e_skin  # rate of energy storage in the skin
        tc_sk = cp_body*alfa*body_weight  # thermal capacity skin
        tc_cr = cp_body*(1-alfa)*body_weight  # thermal capacity core
        d_t_sk = dt*s_skin*body_surface_area/tc_sk  # rate of change skin temperature °C per minute
        d_t_cr = dt*s_core*body_surface_area/tc_cr  # rate of change core temperature °C per minute
        
        # return for solving procedure
        if qsens_return==False:
            return np.asarray([d_t_cr, d_t_sk])
        # return q_sensible
        else:
            return q_sensible
    # dimensionner les vecteurs retour
    Tskin,Tcore=np.zeros(len(time)),np.zeros(len(time))
    Tskin[0],Tcore[0]=t_skin,t_core
    for k,t in enumerate(time):
        # prepare solving with crank-nicolson
        def fc_implicity(Tp, T, theta):
            dTp= fc_TcoreTskin(Tp[0], Tp[1],False,k) # implicit part
            dT = fc_TcoreTskin( T[0], T[1], False,k) # explicit part
            equation = -Tp + T + (1-theta)*dT+ theta*dTp #implicit scheme
            return equation
        
        Tinit=np.asarray([t_core,t_skin])
        Tavant=Tinit
        T_crank=fsolve(fc_implicity, Tinit, args=(Tavant, theta))
        t_core=T_crank[0]
        t_skin=T_crank[1]
        Tskin[k],Tcore[k]=t_skin,t_core
        # get q_sensible
        q_sensible=fc_TcoreTskin(t_core, t_skin, True,k)
        
        t_body = alfa * t_skin + (1 - alfa) * t_core  # mean body temperature, °C
        # sk_sig thermoregulatory control signal from the skin
        sk_sig = t_skin - temp_skin_neutral
        warm_sk = (sk_sig > 0) * sk_sig  # vasodilation signal
        colds = ((-1.0 * sk_sig) > 0) * (-1.0 * sk_sig)  # vasoconstriction signal
        # c_reg_sig thermoregulatory control signal from the skin, °C
        c_reg_sig = t_core - temp_core_neutral
        # c_warm vasodilation signal
        c_warm = (c_reg_sig > 0) * c_reg_sig
        # c_cold vasoconstriction signal
        c_cold = ((-1.0 * c_reg_sig) > 0) * (-1.0 * c_reg_sig)
        # bd_sig thermoregulatory control signal from the body
        bd_sig = t_body - temp_body_neutral
        warm_b = (bd_sig > 0) * bd_sig
        m_bl = (skin_blood_flow_neutral + c_dil * c_warm)/(1 + c_str * colds)
        if m_bl > max_skin_blood_flow:
            m_bl = max_skin_blood_flow
            heat_strain_blood_flow = True
        if m_bl < 0.5: m_bl = 0.5
        m_rsw = c_sw * warm_b * np.exp(warm_sk / 10.7)  # regulatory sweating g/m2/h
        if m_rsw > 500.0: #g/m2/h
            m_rsw = 500.0
            heat_strain_sweating = True
        e_rsw = Lv*m_rsw/3600 # heat lost by vaporization sweat
            
        # calcul e_max
        e_max = (fc_pvs(t_skin) - vapor_pressure)/r_etot
        
        p_rsw = e_rsw / e_max  # ratio heat loss sweating to max heat loss sweating
        w = 0.06 + 0.94 * p_rsw  # skin wetness
        e_diff = w * e_max - e_rsw  # vapor diffusion through skin
        if w > w_max:
            w = w_max
            p_rsw = w_max / 0.94
            e_rsw = p_rsw * e_max
            e_diff = 0.06*(1.0 - p_rsw)*e_max
            heat_strain_w = True
        if e_max < 0:
            e_diff = 0
            e_rsw = 0
            w = w_max
        e_skin = (e_rsw + e_diff)  # total evaporative heat loss sweating and vapor diffusion
        met_shivering = 19.4 * colds * c_cold  # met shivering W/m2
        m = rm + met_shivering
        alfa = 0.0417737 + 0.7451833 / (m_bl + 0.585417)
    q_skin = q_sensible + e_skin  # total heat loss from skin, W
    # p_s_sk saturation vapour pressure of water of the skin
    p_s_sk = fc_pvs(t_skin) 

    # standard environment - where _s at end of the variable names stands for standard
    h_r_s = h_r  # standard environment radiative heat transfer coefficient
    
    h_c_s = 3.0 * np.power(pressure_in_atmospheres, 0.53)
    if not calculate_ce and met > 0.85:
        h_c_met = 5.66 * (met - 0.85) ** 0.39
        h_c_s = max(h_c_s, h_c_met)
    if h_c_s < 3.0:
        h_c_s = 3.0
    
    h_t_s = (h_c_s + h_r_s)  # sum of convective and radiant heat transfer coefficient W/(m2*K)
    r_clo_s = (1.52 / ((met - wme / met_factor) + 0.6944) - 0.1835)  # thermal resistance of clothing, °C M^2 /W
    r_cl_s = 0.155 * r_clo_s  # thermal insulation of the clothing in m2K/W
    f_a_cl_s = 1.0 + k_clo * r_clo_s  # increase in body surface area due to clothing
    f_cl_s = 1.0 / (1.0 + 0.155 * f_a_cl_s * h_t_s * r_clo_s)  # ratio of surface clothed body over nude body
    i_m_s = 0.45  # permeation efficiency of water vapour through the clothing layer
    i_cl_s = (i_m_s * h_c_s / h_t_s * (1 - f_cl_s) / (h_c_s / h_t_s - f_cl_s * i_m_s))  # clothing vapor permeation efficiency
    r_a_s = 1.0 / (f_a_cl_s * h_t_s)  # resistance of air layer to dry heat
    r_ea_s = 1.0 / (lr * f_a_cl_s * h_c_s)
    r_ecl_s = r_cl_s / (lr * i_cl_s)
    h_d_s = 1.0 / (r_a_s + r_cl_s)
    h_e_s = 1.0 / (r_ea_s + r_ecl_s)
    
    def fc_SET(SET):
        return q_skin -h_d_s*(t_skin-SET)-w*h_e_s*(p_s_sk-0.5*fc_pvs(SET))
    _set=fsolve(fc_SET,tdb)[0]
    # print('func',round(_set,2))
    # return Tskin,Tcore,time
    return _set

if __name__ == '__main__':
    
    # Tair, Tmrt,HR, v = 40,25,50, 0.15
    Tair, Tmrt, HR, v = 25,25,50, 0.15
    # Tair, Tmrt,HR, v = 25,25,50, 3
    # Tair, Tmrt,HR, v = 25,25,50, 1.1
    # Tair, Tmrt,HR, v = 25,25,90, 0.15
    T6h=10
    patm= 101325
    met=1.
    iso_9920=True
    iso_9920=False
    SET=set_optimized(Tair,Tmrt,HR,v,T6h,patm,iso_9920)

    print(round(SET,3))

