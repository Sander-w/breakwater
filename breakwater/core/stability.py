import numpy as np

from ..utils.exceptions import NotSupportedError


def xi_critical(Cpl, Cs, P, alpha):
    """ Critical value of the surf-similarity parameter for Van der Meer

    The critical value of the surf-similarity parameter is used in the
    Van der Meer formulas to determine the transition from plunging to
    surging waves.

    Parameters
    ----------
    P : float
        Notional permeability of the structure [-]
    Cpl : float
        Model constant for plunging waves [-]
    Cs : float
        Model constant for surging waves [-]
    alpha : float
        Slope of the structure [rad]

    Returns
    -------
    float
        The critical value of the surf-similarity parameter [-]
    """
    xi_cr = (Cpl/Cs * P**0.31 * np.sqrt(np.tan(alpha)))**(1/(P+0.5))
    return xi_cr


def vandermeer_deep(Hs, Delta, P, Sd, N, xi_m, alpha, safety=1):
    """ Van der Meer formulae for deep water conditions

    For deep water conditions Van der Meer (1988) derived formulae to
    predict the stability of armourstone on uniform armourstone slopes
    with crests above the maximum run-up level. These formulae are only
    valid in deep water conditions, which is defined as
    :math:`h > 3 \\cdot H_{s-toe}`, where :math:`h` is the water depth
    on the toe of the structure and :math:`H_{s-toe}` the significant
    wave height at the toe. The Van der Meer formulae are given by:

    .. math::
       \\frac{H_{s}}{\\Delta D_{n 50}}=c_{p l} P^{0.18}\\left(\\frac{S}
       {\\sqrt{N}}\\right)^{0.2} \\xi_{m}^{-0.5} \\; \\;
       \\text{for plunging waves}

       \\frac{H_{s}}{\\Delta D_{n 50}}=c_{s} P^{-0.13}\\left(\\frac{S}
       {\\sqrt{N}}\\right)^{0.2} \\sqrt{\\cot \\alpha} \\, \\xi_{m}^{P}
       \\; \\; \\text{for surging waves}

    As the formulae are based on a large amount of experiments, there
    is a certain reliability for the model constants. For :math:`C_{pl}`
    this is given by :math:`6.2 \\pm 0.4` and :math:`1 \\pm 0.08` for
    :math:`C_{s}`. Furthermore, each parameter has its own range of
    validity. The range of validity of parameters is presented in the
    table below.

    +------------------------------------+---------------------+-------------+
    | Parameter                          |        Symbol       |  Range      |
    +====================================+=====================+=============+
    | Slope angle                        |      tan(alpha)     | 1:6 - 1:1.5 |
    +------------------------------------+---------------------+-------------+
    | Number of waves                    |          N          |   < 7500    |
    +------------------------------------+---------------------+-------------+
    | Fictitious wave steepness          |        s_om         | 0.01 - 0.06 |
    +------------------------------------+---------------------+-------------+
    | Surf similarity parameter          |        xi_m         |  0.7 - 7    |
    +------------------------------------+---------------------+-------------+
    | Relative buoyant density of armour |        Delta        |    1 - 2.1  |
    +------------------------------------+---------------------+-------------+
    | Relative water depth at toe        |      h/H_s-toe      |    > 3      |
    +------------------------------------+---------------------+-------------+
    | Notional permeability parameter    |          P          |  0.1 - 0.6  |
    +------------------------------------+---------------------+-------------+
    | Armour gradation                   |      Dn85/Dn15      |    < 2.5    |
    +------------------------------------+---------------------+-------------+
    | Damage-storm duration ratio        |      Sd/sqrt(N)     |    < 0.9    |
    +------------------------------------+---------------------+-------------+
    | Stability number                   |   Hs/(Delta Dn50)   |    1 - 4    |
    +------------------------------------+---------------------+-------------+
    | Damage level parameter             |         Sd          |    1 - 20   |
    +------------------------------------+---------------------+-------------+

    Parameters
    ----------
    Hs : float
        The significant wave height, :math:`H_{1/3}` of the incident
        waves at the toe [m]
    Delta : float
        Relative buoyant density [-]
    P : float
        Notional permeability of the structure [-]
    Sd : float
        Damage level parameter [-]
    N : int
        Number of incident waves at the toe of the structure [-]
    xi_m : float
        :math:`\\xi_m`, the surf-similarity parameter computed with
        the mean wave period :math:`T_m` [-]
    alpha : float
        Angle of the front slope [rad]
    safety : float, optional, default: 1
        With this parameter the model constants, Cpl and Cs, can be
        decreased to increase the safety. This is implemented as
        :math:`C_{pl} = \\mu - safety \\cdot \\sigma`.

    Returns
    -------
    Dn50 : float
        the nominal diameter of the armourstone [m]
    """
    Cpl = 6.2 - safety*0.4
    Cs = 1 - safety*0.08

    xi_cr = xi_critical(Cpl, Cs, P, alpha)

    if xi_m < xi_cr:  # Plunging waves
        Dn50 = Hs/(Delta*Cpl*P**0.18
                   * (Sd/np.sqrt(N))**0.2*xi_m**-0.5)
    else:  # Surging waves
        Dn50 = Hs/(Delta*Cs*P**-0.13
                   * (Sd/np.sqrt(N))**0.2*np.sqrt(1/np.tan(alpha))
                   * xi_m**P)
    return Dn50


def vandermeer_shallow(Hs, H2, Delta, P, Sd, N, xi_s_min_1, alpha, safety=1):
    """ Van der Meer formulae for shallow water conditions

    Based on the analysis of the stability of rock-armoured slopes in
    many conditions, mainly focussed on conditions with shallow
    foreshores, it was proposed in Van Gent et al (2003) to modify the
    formulae of Van der Meer (1988) to extend its field of application.
    The Van der Meer formulae for shallow water are:

    .. math::
       \\frac{H_{s}}{\\Delta D_{n 50}}=c_{p l} P^{0.18}\\left(\\frac{S}
       {\\sqrt{N}}\\right)^{0.2}\\left(\\frac{H_{s}}{H_{2 \\%}}\\right)
       \\, \\xi_{m-1,0}^{-0.5} \\; \\; \\text{for plunging waves}

       \\frac{H_{s}}{\\Delta D_{n 50}}=c_{s} P^{-0.13}\\left(\\frac{S}
       {\\sqrt{N}}\\right)^{0.2}\\left(\\frac{H_{s}}{H_{2 \\%}}\\right)
       \\sqrt{\\cot \\alpha} \\, \\xi_{m-1,0}^{P} \\; \\;
       \\text{for surging waves}

    As the formulae are based on a large amount of experiments, there
    is a certain reliability for the model constants. For :math:`C_{pl}`
    this is given by :math:`8.4 \\pm 0.7` and :math:`1.3 \\pm 0.15` for
    :math:`C_{s}`. Furthermore, each parameter has its own range of
    validity. The range of validity of parameters is presented in the
    table below.

    +------------------------------------+---------------------+-------------+
    | Parameter                          |        Symbol       |  Range      |
    +====================================+=====================+=============+
    | Slope angle                        |      tan(alpha)     |  1:4 - 1:2  |
    +------------------------------------+---------------------+-------------+
    | Number of waves                    |          N          |   < 3000    |
    +------------------------------------+---------------------+-------------+
    | Fictitious wave steepness          |        s_om         | 0.01 - 0.06 |
    +------------------------------------+---------------------+-------------+
    | Surf similarity parameter          |        xi_m         |    1 - 5    |
    +------------------------------------+---------------------+-------------+
    | Surf similarity parameter          |      xi_m_min_1     |  1.3 - 6.5  |
    +------------------------------------+---------------------+-------------+
    | Wave height ratio                  |       H2%/Hs        |  1.2 - 1.4  |
    +------------------------------------+---------------------+-------------+
    | Deep-water H over h at the toe     |        Hso/h        | 0.25 - 1.5  |
    +------------------------------------+---------------------+-------------+
    | Armour gradation                   |      Dn85/Dn15      |  1.4 - 2.0  |
    +------------------------------------+---------------------+-------------+
    | Core material - armour ratio       |    Dn50-core/Dn50   |    0 - 0.3  |
    +------------------------------------+---------------------+-------------+
    | Stability number                   |   Hs/(Delta Dn50)   |  0.5 - 4.5  |
    +------------------------------------+---------------------+-------------+
    | Damage level parameter             |         Sd          |      < 30   |
    +------------------------------------+---------------------+-------------+


    Parameters
    ----------
    Hs : float
        The significant wave height, :math:`H_{1/3}` of the incident
        waves at the toe [m]
    H2 : float
        :math:`H_{2\\%}`, wave height exceeded by 2% of the incident
        waves at the toe [m]
    Delta : float
        Relative buoyant density [-]
    P : float
        Notional permeability of the structure [-]
    Sd : float
        Damage level parameter [-]
    N : int
        Number of incident waves at the toe of the structure [-]
    xi_m_min_1 : float
        :math:`\\xi_{s-1.0}`, the surf-similarity parameter computed
        with the energy wave period :math:`T_{m-1.0}` [-]
    alpha : float
        Angle of the front slope [rad]
    safety : float, optional, default: 1
        With this parameter the model constants, Cpl and Cs, can be
        decreased to increase the safety. This is implemented as
        :math:`C_{pl} = \\mu - safety \\cdot \\sigma`.

    Returns
    -------
    Dn50 : float
        the nominal diameter of the armourstone [m]
    """
    Cpl = 8.4 - safety*0.7
    Cs = 1.3 - safety*0.15

    xi_cr = xi_critical(Cpl, Cs, P, alpha)

    if xi_s_min_1 < xi_cr: # Plunging waves
        Dn50 = Hs/(Delta*Cpl* P**0.18
                   * (Sd/np.sqrt(N))**0.2 * (Hs/H2) * xi_s_min_1**-0.5)
    else: # Surging waves
        Dn50 = Hs/(1.6*Cs*P**-0.13*(Sd/np.sqrt(N))**0.2
                   * (Hs/H2) * np.sqrt(1/np.tan(alpha))
                   * xi_s_min_1**P)
    return Dn50


def vandermeer(
        LimitState, Delta, P, N, alpha, slope_foreshore, val='max',
        safety=1, beta = None, logger=None):
    """ Van der Meer formulae for deep and shallow water conditions

    Implementation of Van der Meer formulae for deep and shallow
    water conditions. The function first determines if the water depth
    and wave height are in range of the shallow or deep water formulae,
    and then uses the correct formulae to compute the nominal diameter
    of the armour layer. In case the input is in the range of the deep
    and shallow water formulae both formulae are used to compute the
    diameter, the largest diameter of the two is then returned.

    Parameters
    ----------
    LimitState : :py:class:`LimitState`
        ULS, SLS or another limit state defined with
        :py:class:`LimitState`
    Delta : float
        Relative buoyant density [-]
    P : float
        Notional permeability of the structure [-]
    N : int
        Number of incident waves at the toe of the structure [-]
    alpha : float
        Angle of the front slope [rad]
    slope_foreshore : float
        slope of the foreshore [rad]
    val : {min, max, avg}, optional, default: max
        value to return in case both the deep and shallow water formula
        are valid. min for the lowest value, max for the highest value
        and avg for the average value, default is max.
    safety : float, optional, default: 1
        With this parameter the model constants, Cpl and Cs, can be
        decreased to increase the safety. This is implemented as
        :math:`C_{pl} = \\mu - safety \\cdot \\sigma`.
    logger : dict, optional, default: None
        dict to log messages, must have keys 'INFO' and 'WARNINGS'

    Returns
    -------
    Dn50 : float
        the nominal diameter of the armourstone [m]
    """
    Hs = LimitState.get_Hs('H13')
    H2_per = LimitState.get_H2(slope_foreshore)
    Sd = LimitState['Sd']

    # only the deep water formula is valid
    if LimitState.h/Hs > 3:
        xi_m = LimitState.surf_similarity(alpha=alpha, number='mean')
        Dn50 = vandermeer_deep(
            Hs=Hs, Delta=Delta, P=P, Sd=Sd, N=N, xi_m=xi_m, alpha=alpha,
            safety=safety)
        msg = f'in range of vandermeer_deep with {LimitState.label}'

    # both the deep and shallow water formula are valid
    elif LimitState.h/Hs < 3 and LimitState.h/Hs > 2:
        xi_m = LimitState.surf_similarity(alpha=alpha, number='mean')
        Dn50_d = vandermeer_deep(
            Hs=Hs, Delta=Delta, P=P, Sd=Sd, N=N, xi_m=xi_m, alpha=alpha,
            safety=safety)

        xi_min_1 = LimitState.surf_similarity(alpha=alpha, number='spectral')
        Dn50_s = vandermeer_shallow(
            Hs=Hs, H2=H2_per, Delta=Delta, P=P, Sd=Sd, N=N,
            xi_s_min_1=xi_min_1, alpha=alpha, safety=safety)

        if val == 'max':
            # user wants maximum value
            if Dn50_d >= Dn50_s:
                # deep water formula is normative
                Dn50 = Dn50_d
                formula = 'vandermeer_deep'
            else:
                # shallow water formula is normative
                Dn50 = Dn50_s
                formula = 'vandermeer_shallow'

        elif val == 'min':
            # user wants minimum value
            if Dn50_d >= Dn50_s:
                # shallow water formula is normative
                Dn50 = Dn50_s
                formula = 'vandermeer_shallow'
            else:
                # deep water formula is normative
                Dn50 = Dn50_d
                formula = 'vandermeer_deep'

        elif val == 'avg':
            # compute average
            Dn50 = (Dn50_d + Dn50_s)/2
            formula = 'both'

        else:
            # invalid argument for val
            raise NotSupportedError(
                f'val = {val} is not implemented, must be min, max or avg')


        msg = ('in range of vandermeer_deep and vandermeer_shallow with '
              f'{LimitState.label}, used Dn50 of {formula} as val = {val}')

    # only the shallow water formula is valid
    else:
        xi_min_1 = LimitState.surf_similarity(alpha=alpha, number='spectral')
        Dn50 = vandermeer_shallow(
            Hs=Hs, H2=H2_per, Delta=Delta, P=P, Sd=Sd, N=N,
            xi_s_min_1=xi_min_1, alpha=alpha, safety=safety)
        msg = f'in range of vandermeer_shallow with {LimitState.label}'

    # add msg to log which formula was used
    if logger != None:
        logger['INFO'].append(msg)


    return Dn50


def hudson(H, Kd, Delta, alpha):
    """ Hudson formula for armour layer stability (Hudson, 1959)

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
    Kd : int
        Stability coefficient [-]
    Delta : float
        Relative buoyant density [-]
    alpha : float
        Angle of the front slope [rad]

    Returns
    -------
    Dn50 : float
        the nominal diameter [m]
    """
    Dn50 = H/(Delta*(Kd/np.tan(alpha))**(1/3))
    return Dn50

def vangent(Dn50_per, beta, cb):
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
    Dn50_per: float
        Dn50 as calculated for stability under perpendicular wave attack [m]
    beta: float
        Angle of oblique wave attack [rad]
    cb: float
        Correction factor for wave type [-]
        cb = 0.42 for rock slopes with short-crested waves
        cb = 0.35 for rock slopes with long-crested waves
        cb = 0.35 for cubes in a double layer
        cb = 0    for cubes in a single layer
    
    """
    
    Dn50_beta = Dn50_per*((1-cb)*np.cos(beta)**2+cb)
    
    return(Dn50_beta)