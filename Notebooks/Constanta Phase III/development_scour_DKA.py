# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 09:42:11 2022

@author: JY6
"""

import numpy as np
from breakwater.utils.exceptions import user_warning

#%% Input for testing

g = 9.81
Hs = 2.24
Tp = 8.14
h = 7.1
C2 = 1

#%% Actual function
def sumer_fredsoe(Hs, Tp, g, h, C2 = 1):
    """ Compute the equilibrium scour depth S with the empirical expression
    of Sumer and Fredsøe (2002)

        .. math::
           S = 0.01 * C_{2} * \left( 
               \frac{T_{p}*\sqrt{g*H_{s}}}{h} 
               \right)^{1.5}*H_{s}

        Computes the equilibrium scour depth for a given peak period and
        significant wave height.

        Parameters
        ----------
        Hs : float
            significant wave height [m]
        Tp   : float
            peak wave period [s]
        g   : float
            gravity constant [m/s^2]
        h   : float
            water depth [m]
        C2  : float
            Uncertainty factor, 1 by default
        
        Returns
        -------
        S : float
            the equilibrium scour depth [m]

        Raises
        ------
        InputError
            If one of the inputs is missing
        """
        
    S = 0.01*C2*((Tp*np.sqrt(g*Hs))/(h))**1.5*Hs


    return S

def v_scour(S, t, cot_a_scour_hole):
    """ Compute added toe volume necessary for mitigating scour effects.

        .. math::
           V  = t * \sqrt{S^2+(S*cot(\alpha))^2}

        Compute the extra necessary toe volume for mitigating the effects of
        scour. The formula assumes that the extra volume will act as a 
        falling apron.

        Parameters
        ----------
        S    : float
            Scour depth [m]
        t    : float
            Layer thickness after falling [Dn50]
        cota : float
            Scour hole slope [-]
        
        Returns
        -------
        V : float
            Toe volume necessary to mitigate the effects of 
            the input scour [m^3/m]

        Raises
        ------
        InputError
            If one of the inputs is missing
        """
        
    V = t * np.sqrt(S**2+(S*cot_a_scour_hole)**2)
        
    return V

#%% Test calculations
S = sumer_fredsoe(Hs, Tp, g, h, C2)