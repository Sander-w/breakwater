# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:38:44 2022

@author: JY6
"""

import numpy as np
from breakwater.utils.exceptions import user_warning

def surf_similarity(tana, H, T, g):
        """ Compute the surf similarity parameter

        .. math::
           \\xi = \\frac{tan{\\alpha}}{\\sqrt{\\frac{2*pi*H}{g*T^2}}}

        Computes the surf similarity parameter, also known as the
        Iribarren number. The correct wave height and period need to be input
        in the equation

        Parameters
        ----------
        tana : float
            tangent slope of the structure [rad]
        H   : float
            wave height, usually Hm0 or Hs
        T   : float
            wave period, usually Tm-1,0
        g   : float
            gravity constant
        
        Returns
        -------
        xi : float
            the surf similarity parameter [-]

        Raises
        ------
        InputError
            If one of the inputs is missing
        """
        

        xi = tana/np.sqrt((2*np.pi/g)*(H/T**2))


        return xi

def calc_beta(Dir_structure, Dir_wave):
    """ Compute the wave obliqueness angle
    
    Computes the angle of wave obliqueness of a cross-section.
    Directions are given in the nautical convention.

    Parameters
    ----------
    Dir_structure : float
        Direction of the structure normal in nautical convention. 0 deg is
        facing north, 180 deg is facing south. Negative input is also 
        accepted.
        [deg]
    Dir_wave : float
        Direction where the waves come from, in nautical convention. 0 deg is
        coming from north, 90 from east, 180 from south. Negative input is
        also accepted.
        [deg]

    Returns
    -------
    beta : float
        Angle between structure normal and wave direction [deg]

    """
    
    # Calculate beta. Calculation method based on 
    # AWP-#3504459-V48-IHS-COA-xxx-CAL_Armour_Stability_under_Waves.XLSM
    if (Dir_wave-Dir_structure)%360 < 180:
        beta = (Dir_wave-Dir_structure)%360
    else: 
        beta = (Dir_structure-Dir_wave)%360
    
    
    
    # Give warning if input is larger than 360
    if abs(Dir_wave)>360:
        user_warning(f'Wave direction is {Dir_wave} degrees. '
                     'Make sure this is not an input mistake.')
    elif abs(Dir_structure)>360:
        user_warning(f'Structure direction is {Dir_structure} degrees. '
                     'Make sure this is not an input mistake.')
    
    
    
    return(beta)



def eurotop2018_6_5(
        g, Hm0, q, gamma_f, gamma_beta,
        safety=0,
        Gc=None, Dn50=None, limit=True):
    """ Compute the crest freeboard of a rubble mound breakwater

    Computes the crest freeboard of a rubble mound breakwater using
    equation 6.5 from EurOtop (2018).

    .. math::
       \\frac{q}{\\sqrt{g \\cdot H_{m 0}^{3}}}=0.09 \\cdot \\exp \\left
       [-\\left(1.5 \\frac{R_{c}}{H_{m 0} \\cdot \\gamma_{f} \\cdot
       \\gamma_{\\beta} \\cdot \\gamma^{*}}\\right)^{1.3}\\right]

    The reliability of the first equation is given by
    :math:`\\sigma` (0.09) = 0.0135 and
    :math:`\\sigma` (1.5) = 0.15.

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    xi_m_min_1 : float
        :math:`\\xi_{m-1.0}`, the surf-similarity parameter computed
        with the energy wave period :math:`T_{m-1.0}` [-]
    alpha : float
        Angle of the front slope [rad]
    gamma_f : float
        influence factor for roughness [-]
    gamma_beta : float
        influence factor for obliqueness [-]
    safety : float, optional, default: 0
        safety factor of the design, positive values increase the safety
        of the design by increasing the mean value of the model constants
        with the number of standard deviations specified. In accordance
        with the recommendation from EurOtop (2018), the default value is
        set to 1 standard deviation.
    Gc : float, optional, default: None
        Effect of armoured crest width [m]
    Dn50 : float, optional, default: None
        nominal diameter of the armour units on the crest [m]
    limit : bool, optional, default: True
        If True, the discharge will be set to the limit for zero
        discharge in case the given discharge is below this limit. If
        False, the discharge will not be changed.

    Returns
    -------
    Rc : float
        the crest freeboard of the structure [m]
    """
    
    gam_f = gamma_f
    gam_beta = gamma_beta

    constant_c = 0.09 + safety*0.0135
    constant_d = 1.5 - safety*0.15

    # check given overtop discharge against EurOtop limit for zero discharge
    limit_zero_q = np.sqrt(g*Hm0**3) * 0.01
    if limit_zero_q > q:
        msg = (f'Given value of q = {q}, is below EurOtop limit for zero '
               f'discharge, which is {np.round(limit_zero_q, 3)}')

        if limit:
            msg += f'. Using {np.round(limit_zero_q, 3)} instead as limit=True'
            q = limit_zero_q

        user_warning(msg)

    # in case of a wide crest the overtopping in the equation can be increased
    if Gc != None and Gc > 3*Dn50:
        Cr = 3.06 * np.exp(-1.5*Gc/Hm0)
        q = q / min(Cr, 1)

    q = q/1000

    # compute arguments of the ln (to check if arg < 1)
    arg = q/(constant_c*np.sqrt(g*Hm0**3))



    if arg <= 1: # handle for RuntimeWarning
        Rc = ((Hm0*gam_f*gam_beta/constant_d)
                    * (-np.log(arg))**(1/1.3))
    #    Rc = (((-1)*(np.log(q/(np.sqrt(g*Hm0**3)*0.09))))**(1/1.3)*Hm0*gamma_f*gamma_beta/1.5)
    else:
        Rc = 0
        user_warning('Encountered negative freeboard')



    return Rc

def gamma_beta_eurotop_2018_6_9(beta):
    """ 
    Adaptation to be in line with Constanta Phase II. Based on 
    breakwater.core.overtopping.gamma_beta
    
    Influence factor for oblique wave attack

    Computes the influence factor for oblique wave attack with equation
    5.29 from EurOtop (2018).

    .. math::
       \\gamma_{\\beta} = 1 - 0.0063 \\mid \\beta \\mid

    with a maximum of :math:`\\gamma_{\\beta} = 0.496` for
    :math:`\\mid \\beta \\mid > 80`
    
    if :math: `\\beta > 0.496`, :math:`\\gamma_{\\beta} = 0.01`

    Parameters
    ----------
    beta : float
        the angle of wave attack [deg]

    Returns
    -------
    gamma_beta : float
        the influence factor for oblique wave attack
    """
    #beta_deg = beta*180/np.pi
    beta_deg = beta

    if beta_deg <= 80:
        return 1 - 0.0063 * np.abs(beta_deg)
    elif beta_deg > 90:
        return 0.001
    else:
        return 0.496
    
def eurotop_6_5_q(
        g, Hm0, Rc, gamma_f, gamma_beta,
        safety=0,
        Gc=None, Dn50=None, limit=True):
    """ Compute the overtopping over a rubble mound breakwater

    Computes the average overtopping discharge of a rubble mound breakwater 
    using equation 6.5 from EurOtop (2018).

    .. math::
       q = \sqrt{g*H_{m0}^3}*0.09*exp[-(1.5\frac{R_c}{H_{m0}*
                      \gamma_f*\gamma_\beta})^{1.3}]

    The reliability of the first equation is given by
    :math:`\\sigma` (0.09) = 0.0135 and
    :math:`\\sigma` (1.5) = 0.15.

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    Rc : float
        freeboard [m]
    xi_m_min_1 : float
        :math:`\\xi_{m-1.0}`, the surf-similarity parameter computed
        with the energy wave period :math:`T_{m-1.0}` [-]
    alpha : float
        Angle of the front slope [rad]
    gamma_f : float
        influence factor for roughness [-]
    gamma_beta : float
        influence factor for obliqueness [-]
    safety : float, optional, default: 0
        safety factor of the design, positive values increase the safety
        of the design by increasing the mean value of the model constants
        with the number of standard deviations specified. In accordance
        with the recommendation from EurOtop (2018), the default value is
        set to 1 standard deviation.
    Gc : float, optional, default: None
        Effect of armoured crest width [m]
    Dn50 : float, optional, default: None
        nominal diameter of the armour units on the crest [m]
    limit : bool, optional, default: True
        If True, the discharge will be set to the limit for zero
        discharge in case the given discharge is below this limit. If
        False, the discharge will not be changed.

    Returns
    -------
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    """
    
    gam_f = gamma_f
    gam_beta = gamma_beta

    constant_c = 0.09 + safety*0.0135
    constant_d = 1.5 - safety*0.15

    # compute arguments of the ln (to check if arg < 1)
    # arg = q/(constant_c*np.sqrt(g*Hm0**3))

    # if arg <= 1: # handle for RuntimeWarning
    #     Rc = ((Hm0*gam_f*gam_beta/constant_d)
    #                 * (-np.log(arg))**(1/1.3))
    # #    Rc = (((-1)*(np.log(q/(np.sqrt(g*Hm0**3)*0.09))))**(1/1.3)*Hm0*gamma_f*gamma_beta/1.5)
    # else:
    #     Rc = 0
    #     user_warning('Encountered negative freeboard')

    q = (np.sqrt(g*Hm0**3)
         *constant_c
         *np.exp(-(constant_d*(Rc)/(Hm0*gam_f*gam_beta))**1.3)
         )
    
    q = q*1000

    # check given overtop discharge against EurOtop limit for zero discharge
    limit_zero_q = np.sqrt(g*Hm0**3) * 0.01
    if limit_zero_q > q:
        msg = (f'Given value of q = {q}, is below EurOtop limit for zero '
               f'discharge, which is {np.round(limit_zero_q, 3)}')

        if limit:
            msg += f'. Using {np.round(limit_zero_q, 3)} instead as limit=True'
            q = limit_zero_q

        user_warning(msg)
        
        # in case of a wide crest the overtopping in the equation can be increased
    if Gc != None and Gc > 3*Dn50:
        Cr = 3.06 * np.exp(-1.5*Gc/Hm0)
        q = q * min(Cr, 1)

    

    return q