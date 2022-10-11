# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 08:25:42 2022

@author: JY6
"""

import numpy as np
from breakwater.utils.exceptions import user_warning, NotSupportedError

def vangent_armour_reduction(beta, cb, beta_max):
    """ Van Gent formula for reduction of stone diameter under 
    oblique wave attack (Van Gent, 2014)
    
    The Van Gent formula is based on experimental research on the
    stability of rock on a slope under oblique wave impact. The 
    formula is derived by fitting a line through data points. 
    The formula calculates a reduced Dn50 for oblique attack.
    
    .. math::
        Dn_{50 \beta} = Dn50_per{(1-cb)*cos(beta)^2+cb}
    
    .. warning::
        This function is under development. The validity limits have not
        been established or added.
        
    .. note:
        
        
    Parameters
    ---------
    beta: float
        Angle of oblique wave attack [deg]
    cb: float
        Correction factor for wave type [-]
        cb = 0.42 for rock slopes with short-crested waves
        cb = 0.35 for rock slopes with long-crested waves
        cb = 0.35 for cubes in a double layer
        cb = 0    for cubes in a single layer
    beta_max: float
        Angle above which no armour is required [deg]
    
    """
    if beta <= 90:
        beta = beta/180*(np.pi)
        gamma_beta = (1-cb)*np.cos(beta)**2+cb
    elif beta <= beta_max:
        gamma_beta = cb
    else:
        gamma_beta = 0.001
    
    return gamma_beta

def hudson_fixed_slope(H, Dv, Delta):
    """ Hudson formula for armour layer stability (Hudson, 1959), 
    applied for a fixed slope

    The Hudson formula is based on experimental research on the stability
    of rock and armour units (Tetrapods and Dolose) on a slope. The
    formula is derived by fitting a line through the data points. All
    unknown factors are included in the 'dustbin' factor :math:`k_D`
    (Van den Bos and Verhagen, 2018). The Hudson formula expressed in
    the form with the stability number is:

    .. math::
       \\frac{H}{\\Delta D_{n 50}}=\\sqrt[3]{k_{D} \\cot \\alpha}

    .. warning::
       The formula is based on regular waves (Hudson, 1959), therefore
       the wave height H to use in case of irregular waves is not
       defined. CERC (1984, p. 7-203) advised in the Shore Protection
       Manual to use :math:`H_{1/10}`, the average of the highest 10%
       of the waves.

    .. note::
       Although the Hudson formula is derived for the stability of rock
       and armour units, Van den Bos and Verhagen (2018) recommend to
       use the Van der Meer formula for the stability of rock. As the
       Van der Meer formula includes the effects of storm duration,
       wave period, the structures permeability and a clearly defined
       damage level (CIRIA, CUR, CETMEF, 2007, p. 567).

    Parameters
    ----------
    H : float
        The design wave height [m]
    Dv : int
        Design value [-] based on tability coefficient [-] and slope [-]
    Delta : float
        Relative buoyant density [-]
    

    Returns
    -------
    Dn50 : float
        the nominal diameter [m]
    """
    Dn50 = H/(Delta*Dv)
    
    return Dn50

def gamma_f_vangent_pozueta(roughness):
    """
    Determination of roughness coefficient for rear side armour stability,
    based on slope classification for equations 5.195, 5.196 and 5.197 in the
    Rock Manual
    """
    
    #Find value for gamma_f based on roughness input
    if roughness == 'rough':
        gamma_f = 0.55
    elif roughness == 'smooth':
        gamma_f = 1
    else:
        raise NotSupportedError(
            (f"Unknown input for roughness: {roughness}. "
             "Input either 'rough' or 'smooth' ")) 

    return gamma_f


def rock_manual_5_194(Sd, 
                      N, 
                      u_1_percent, 
                      Tm_min_1,
                      Delta,
                      cota_rear,
                      Rc_rear,
                      Hs):
    
    """ Van Gent Pozueta (2005) method for rear side armour stability.
    
    The Van Gent and Pozueta formula is based on experimental research on
    the stability of rock on the rear side of a breakwater. The formula is 
    derived by fitting a line through data points. The formula calculates a 
    reduced Dn50 for rear side armour.
    
    .. math::
        D_{n 50} = 0.008 \left( \frac{S_d}{\sqrt{N}} \right)^{-1/6}
        \left( \frac{u_{1\%}T_{m-1,0}}{\sqrt{\Delta}}\right)
        \left(cot \alpha_{rear} \right)^{-2.5/6}
        \left(1+10exp \left( -R_{c,rear}/H_s \right) \right)^{1/6}

    .. warning::
       This function is under development. The validity limits have not
        been added. The function has not been tested for other cases than
        Constanta Phase III

    .. note::
       
        
    Parameters
    ----------
    Sd : float
        Allowed damage number [-]
    N : int
        Number of waves in the design storm [-]
    Hs : float
        Significant wave height [m]
    T_m_min_1 : float
        Energy wave period [s]
    cota_rear: float
        Rear side slope [-]
    Rc_rear : float
        Crest freeboard relative to water level [m]
    Delta : float
        Relative buoyant density [-]
    u_1_percent : float
        
    

    Returns
    -------
    Dn50 : float
        the nominal required diameter of rear side armour stone [m]
    
    """
    
    Dn50 = (0.008 
            * (Sd/np.sqrt(N))**(-1/6) 
            * ((u_1_percent*Tm_min_1)/np.sqrt(Delta)) 
            *  (cota_rear)**(-2.5/6) 
            * (1+10*np.exp(-Rc_rear/Hs))**(1/6))
    
    
    
    return Dn50

# Sd = 2 
# N = 3407
# u_1_percent = 5.7722
# Tm_min_1 = 7.39
# Delta = 1.574
# cota_rear =1.5
# Rc_rear = 2.18
# Hs = 2.67

# Dn50_rear = rock_manual_5_194(Sd, N, u_1_percent, Tm_min_1, Delta, cota_rear, Rc_rear, Hs)

def rock_manual_5_195(g,
                      slope_roughness,
                      crest_roughness,
                      Ru_1_percent,
                      Rc,
                      B,
                      Hs):
    """ calculation of u1% as input for Van Gent and Pozueta rear side
    armour stability 'rock_manual_5_194'
    
    .. math::
        u_{1\%} = 1.7 \left. \left(g \gamma_{f-c} \right)^{0.5}
        \left( \frac{R_{u1\%}-R_c}{\gamma_{f}} \right)^{0.5}
        \middle/ \left( 1+0.1\frac{B}{H_s} \right) \right.

    .. warning::
       This function is under development. The validity limits have not
        been added. The function has not been tested for other cases than
        Constanta Phase III

    .. note::
       
        
    Parameters
    ----------
    g : float
        Gravity constant [m/s2]
    slope_roughness : str
        Indication of seaward slope roughness based on the classification 
        in the Rock Manual.
        'rough' for armourstone slopes, 
        'smooth' for impermeable slopes
    crest_roughness : str
        Indication of crest roughness based on the classification in the
        Rock Manual.
        'rough' for armourstone crests
        'smooth' for impermeable crests.
    Ru_1_percent : float
        Run-up level exceeded by 1 percent of incident waves [m]
    Rc : float
        Crest height [m]
    B : float
        Crst width [m]
    Hs : float
        Significant wave height [m]

    Returns
    -------
   u_1_percent : float
        maximum velocity (depth-averaged) at the rear side of the crest [m/s]
    
    """
    gamma_f_crest = gamma_f_vangent_pozueta(slope_roughness)
    gamma_f_for_rear = gamma_f_vangent_pozueta(slope_roughness)
    
    u_1_percent = ((1.7
                   *(g*gamma_f_crest)**0.5
                   *((Ru_1_percent-Rc)/gamma_f_for_rear)**0.5)
                   /
                   (1 +0.1*B/Hs))
      
    return u_1_percent



# g = 9.81
# slope_roughness = 'rough'
# crest_roughness = 'rough'
# Ru_1_percent = 4.05
# Rc = 2.18
# B = 7.01
# Hs = 2.67


# u_1_percent = rock_manual_5_195(g, 
#                   slope_roughness, 
#                   crest_roughness, 
#                   Ru_1_percent, 
#                   Rc, B, Hs)

# print(f'u1%: {u_1_percent:.2f}')

def rock_manual_5_196(Hs, xi_s_min_1, slope_roughness, beta):
    """ Calculation of Ru1% as input for Van Gent and Pozueta rear side
    armour stability 'rock_manual_5_195'. The function chooses wether to use
    equation 5.196 or 5.197 automatically.
    
    .. math::
        u_{1\%} = 1.7 \left. \left(g \gamma_{f-c} \right)^{0.5}
        \left( \frac{R_{u1\%}-R_c}{\gamma_{f}} \right)^{0.5}
        \middle/ \left( 1+0.1\frac{B}{H_s} \right) \right.

    .. warning::
       This function is under development. The validity limits have not
        been added. The function has not been tested for other cases than
        Constanta Phase III

    .. note::
       
        
    Parameters
    ----------
    Hs : float
        Significant wave height [m]
    xi_s_min_1 : float
        Surf similarity parameter based on Hs and Tm_min_1
    slope_roughness : str
        Indication of rougness of seaward slope for determination 
        of gamma_f. 
        'rough' for armourstone slopes 
        'smooth' for impermable slopes
    beta : float
        Angle of wave incidence used for calculation of gamma_beta [deg]

    Returns
    -------
    gamma_f : float
        Roughness reduction coefficient [-]
    gamma_beta : float
        Obliqueness reduction coefficient [-]
    gamma : float
        Combined reduction coefficient [-]
    Ru_1_percent : float
        maximum velocity (depth-averaged) at the rear side of the crest [m/s]
    
    """
    
    
    
    # Setting of local input coefficients, see Rock Manual p 632
    c0 = 1.45
    c1 = 5.1
    c2 = 0.25*c1**2/c0
    p = 0.5*c1/c0
    
    # Slope roughness coefficient based on roughness classification
    gamma_f = gamma_f_vangent_pozueta(slope_roughness)
    
    # Obliqueness coefficient based on oblique wave angle
    if beta <= 80:
        gamma_beta = 1-0.0022*beta
    elif beta > 80:
        gamma_beta = 1-0.0022*80
    else: 
        user_warning("Unexpected beta input")
    
    # Total reduction coefficient
    gamma = gamma_f*gamma_beta
    
    # Calculate Ru1% according to 5.196 or 5.197
    if xi_s_min_1 <= p:
        Ru_1_percent = (c0*xi_s_min_1) * (gamma*Hs)
    else:
        Ru_1_percent = (c1-c2/xi_s_min_1) * (gamma*Hs)
    
    
    # # Print intermediate values for development purposes
    # print(f'p         : {p:.2f}')
    # print(f'c2        : {c2:.2f}')
    # print(f'gamma     : {gamma:.2f}')
    # print(f'gamma_beta: {gamma_beta:.2f}')
    # print(f'gamma_f   : {gamma_f:.2f}')
    
    
    return Ru_1_percent, gamma_f, gamma_beta, gamma


# Hs = 2.67
# xi_s_min_1 = 3.7674
# slope_roughness = 'rough'
# beta = 13.3

# Ru_1_percent = rock_manual_5_196(Hs, xi_s_min_1, slope_roughness, beta)

# print(f'Ru1% : {Ru_1_percent:.2f}')