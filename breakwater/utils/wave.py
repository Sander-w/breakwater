import numpy as np
from scipy.optimize import fsolve

def shoaling_coefficient(h, T, H0, linear=False):
    """ Closed form equation for the non-linear shoaling coefficient

    Computes the non-linear shoaling coefficient using a closed form
    equation. The equation uses the coefficient derived in Kweon and
    Goda (1996).

    Parameters
    ----------
    h : float
        water depth [m]
    T : float
        wave period, Goda (2000) advises to use :math:`T_{1/3}` [s]
    H0 : float
        equivalent deep water significant wave height [m]
    linear : bool, optional, default: False
        if true the linear shoaling coefficient is returned, if false
        the non-linear shoaling coefficient is returned

    Returns
    -------
    Ks : float
        the non-linear shoaling coefficient
    """
    # set gravitational acceleration
    g = 9.81

    # compute wave length with the dispersion relation
    L = dispersion(T=T, h=h)
    L0  = g*T**2/(2*np.pi)

    # compute wave number
    k = 2*np.pi/L

    # compute ratio c_g (group velocity) to c (wave celerity)
    n = 0.5*(1 + 2*k*h/np.sinh(2*k*h))

    # compute the linear shoaling coefficient
    Ksi = np.sqrt(1/np.tanh(k*h) * 1/(2*n))

    if linear:
        return Ksi
    else:
        # compute the non-linear shoaling coefficient
        K = Ksi + 0.0015 * (h/L0)**-2.87 * (H0/L0)**1.27
        return K


def dispersion(T, h):
    """ Dispersion relationship

    uses fsolve to find the wave length

    Parameters
    ----------
    T : float
        wave period [s]
    h : float
        water depth [m]

    Returns
    L : float
        the wave length [m]
    """
    g = 9.81

    # deep water wave length
    L0  = g*T**2/(2*np.pi)

    # define dispersion relation with a lambda function to solve L
    dispersion = lambda L : L0 * np.tanh(2*np.pi*h/L) - L

    # solve for L
    L = fsolve(dispersion, L0)[0]

    return L
