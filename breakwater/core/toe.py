import numpy as np

from ..utils.wave import dispersion

def toe_stability(Hs, h, ht, Delta, Nod):
    """ Toe stability formula of Van der Meer (1998)

    The armour layer should is supported by a toe, the formula of Van
    der Meer (1998) computes the minimal required nominal diameter of
    the stones necessary. The formula is given as:

    .. math::
       \\frac{H_{s}}{\\Delta D_{n50}}=\\left(6.2 \\frac{h_{t}}{h}+2
       \\right) N_{o d}^{0.15}

    The formula is based on experiments and the range of validity of
    the parameters can be seen in the table below.

    +---------------------------------+---------+-----------+
    | Parameter                       | Symbol  | Range     |
    +=================================+=========+===========+
    | toe depth over water depth      | ht/h    | 0.4 - 0.9 |
    +---------------------------------+---------+-----------+
    | toe depth over nominal diameter | ht/Dn50 |   3 - 25  |
    +---------------------------------+---------+-----------+

    Parameters
    ----------
    Hs : float
        The significant wave height of the incident waves [m]
    h : float
        The water depth in front of the toe [m]
    ht : float
        The water depth on top of the toe [m]
    Delta : float
        Relative buoyant density [-]
    Nod : float
        Damage number [-]

    Returns
    -------
    Dn50 : float
        Nominal diameter of the armourstones in the toe [m]
    """
    Dn50 = Hs/(Delta*(2 + 6.2*(ht/h)**2.7)*Nod**0.15)
    return Dn50


def toe_berm_stability(Hs, T, d, Bm, Delta, beta=0, alpha_s=0.45):
    """ Berm protection to caisson or vertical wall breakwaters

    Compute the nominal diameter of the armourstone on the berm of a
    caisson or vertical wall breakwaters with the modified Tanimoto
    formula from Takahasi (2002).

    .. math::
       D_{n50} = \\frac{H_{s}}{\\Delta N_{s}}

    in which:

    .. math::
       N_{s}=\\max \\left\\{1.8,\\left(1.3 \\frac{1-\\kappa}{\\kappa^{1 / 3}}
       \\frac{h^{\\prime}}{H_{s}}+1.8 \\exp \\left[-1.5 \\frac{
       (1-\\kappa)^{2}}{\\kappa^{1 / 3}} \\frac{h^{\\prime}}{H_{s}}
       \\right]\\right)\\right\\}

    where:

    .. math::
       \\kappa=\\frac{4 \\pi h^{\\prime} / L^{\\prime}}{\\sinh \\left(4
       \\pi h^{\\prime} / L^{\\prime}\\right)} \\kappa_{2}

    and:

    .. math::
       \\kappa_{2}=\\max \\left\{\\alpha_{S} \\sin ^{2} \\beta \\cos ^{2}
       \\left(\\frac{2 \\pi x}{L^{\\prime}} \\cos \\beta\\right)\\right.,
       \\left.\\cos ^{2} \\beta \\sin ^{2}\\left(\\frac{2 \\pi x}
       {L^{\\prime}} \\cos \\beta\\right)\\right\\}


    Parameters
    ----------
    Hs : float
        mean of the highest 1/3 of the wave heights [m]
    T : float
        wave period (s)
    d : float
        water depth at which the armour is placed [m]
    Bm : float
        width of the berm [m]
    Delta : float
        Relative buoyant density [-]
    beta : float, optional, default: 0
        angle between direction of wave approach and a line normal to
        the breakwater [rad]
    alpha_s : float, optional, default: 0.45
        factor the include the effect of the slope, the value of 0.45 is
        given by measures data (Goda, 2000)

    Returns
    -------
    Dn50 : float
        Nominal diameter of the armourstones on the berm [m]
    """
    L = dispersion(T=T, h=d)
    kh = 2*np.pi*d/L

    kappa2 = max([alpha_s*np.sin(beta)**2 * np.cos(2*np.pi*Bm*np.cos(beta)/L)**2,
                  np.cos(beta)**2 * np.sin(2*np.pi*Bm*np.cos(beta)/L)**2])

    kappa = 2*kh/(np.sinh(2*kh)) * kappa2

    a = (1-kappa)/kappa**(1/3)
    Ns = (1.3 * a * d/Hs + 1.8 * np.exp(-1.5*d*a*(1-kappa)/Hs))
    Ns = max([1.8, Ns])

    Dn50 = Hs/(Delta*Ns)

    return Dn50
