# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 08:25:42 2022

@author: JY6
"""

import numpy as np

def vangent_armour_reduction(beta, cb, beta_max):
    """ Van Gent formula for reduction of stone diameter under 
    oblique wave attack (Van Gent, 2014)
    
    The Van Gent formula is based on experimental research on the
    stability of rock on a slope under oblique wave impact. The 
    formula is derived by fitting a line through data points. 
    The formula calculates a reduced Dn50 for oblique attack.
    
    .. math::
        Dn50_beta = Dn50_per{(1-cb)*cos(beta)^2+cb}
    
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
        beta = beta/360*(2*np.pi)
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