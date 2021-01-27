import numpy as np

from ..utils.exceptions import user_warning, NotSupportedError

def gamma_f(
        armour_layer, xi_m_min_1, layers=None, permeability=None,
        placement=None):
    """ Influence factor for the permeability and roughness of the slope

    Computes the influence factor on roughness with table 6.2 from
    EurOtop (2018). These values are derived for 2.8
    :math:`\\leq \\xi_{m-1.0} \\leq` 4.5, in case of larger surf-similarity
    parameters the influence factor for roughness is increased using
    equation 6.7 from EurOtop (2018)

    +--------------------------------------+-------------------+
    | Type of armour layer                 | :math:`\\gamma_f`  |
    +======================================+===================+
    | Smooth impermeable surface           |       1.00        |
    +--------------------------------------+-------------------+
    | Rock (1 layer, impermeable core)     |       0.60        |
    +--------------------------------------+-------------------+
    | Rock (1 layer, permeable core)       |       0.45        |
    +--------------------------------------+-------------------+
    | Rock (2 layers, impermeable core)    |       0.55        |
    +--------------------------------------+-------------------+
    | Rock (2 layers, permeable core)      |       0.40        |
    +--------------------------------------+-------------------+
    | Cubes (1 layer, flat positioning)    |       0.49        |
    +--------------------------------------+-------------------+
    | Cubes (2 layers, random positioning) |       0.47        |
    +--------------------------------------+-------------------+
    | Antifers                             |       0.50        |
    +--------------------------------------+-------------------+
    | HARO                                 |       0.47        |
    +--------------------------------------+-------------------+
    | Tetrapods                            |       0.38        |
    +--------------------------------------+-------------------+
    | Dolos                                |       0.43        |
    +--------------------------------------+-------------------+
    | Accropode I                          |       0.46        |
    +--------------------------------------+-------------------+
    | Xbloc, CoreLoc, Accropode II         |       0.44        |
    +--------------------------------------+-------------------+
    | XblocPlus                            |       0.45        |
    +--------------------------------------+-------------------+
    | Cubipods (1 layer)                   |       0.49        |
    +--------------------------------------+-------------------+
    | Cubipods (2 layers)                  |       0.47        |
    +--------------------------------------+-------------------+


    Parameters
    ----------
    armour_layer : str
        name of the material of the armour layer, supported materials:
        Smooth, Rock, Cubes, Antifers, HARO, Tetrapods, Dolos,
        Accropode I, Xbloc, XblocPlus, CoreLoc, Accropode II, Cubipods
    xi_m_min_1 : float
        the surf-similarity parameter [-]
    layers : {1, 2}, optional, default: None
        number of layers in the armour layer, required if the armour
        layer is made out of: Rock, Cubes or Cubipods
    permeability : {'permeable', 'impermeable'}, optional, default: None
        permeability of the armour layer, required if the armour layer
        is made out of Rock
    placement : {'flat', 'random'}, optional, default: None
        placement of the armour layer, required if the armour layer is
        made out of Cubes

    Returns
    -------
    gamma_f : float
        the influence factor for roughness

    Raises
    ------
    Keyerror
        If the armour layer is not in table 6.2 from EurOtop (2018)
    """
    table = {'Smooth': 1,
             'Rock - 1 - impermeable': 0.6,
             'Rock - 1 - permeable': 0.45,
             'Rock - 2 - impermeable': 0.55,
             'Rock - 2 - permeable': 0.4,
             'Cubes - 1 - flat': 0.49,
             'Cubes - 2 - random': 0.47,
             'Antifers': 0.5,
             'HARO': 0.47,
             'Tetrapods': 0.38,
             'Dolos': 0.43,
             'Accropode I': 0.46,
             'Xbloc': 0.44,
             'XblocPlus': 0.45,
             'CoreLoc': 0.44,
             'Accropode II': 0.44,
             'Cubipods - 1': 0.49,
             'Cubipods - 2': 0.47}
    if armour_layer == 'Smooth':
        key = 'Smooth'
    elif armour_layer == 'Rock':
        key = f'Rock - {layers} - {permeability}'
    elif armour_layer == 'Cubipods':
        key = f'Cubipods - {layers}'
    elif armour_layer == 'Cubes':
        key = f'Cubes - {layers} - {placement}'
    else:
        key = f'{armour_layer}'

    # check if the armour layer (key) is in the table
    if key in table.keys():
        # get the gamma_f value
        gamma_f = table[key]
    else:
        # raise error as the given armour layer is not supported
        raise NotSupportedError(
            ('EurOtop (2018) does not support a roughness factor for '
             f'{armour_layer}'))

    if xi_m_min_1 > 5:
        gamma_f = gamma_f + (xi_m_min_1-5)*(1-gamma_f)/5
    elif xi_m_min_1 < 2.8:
        user_warning(
            (f'xi out of range in gamma_f, {np.round(xi_m_min_1,2)} < 2.8. '
             f'Continued with gamma_f = {gamma_f}'))
    return gamma_f

def gamma_beta(beta):
    """ Influence factor for oblique wave attack

    Computes the influence factor for oblique wave attack with equation
    5.29 from EurOtop (2018).

    .. math::
       \\gamma_{\\beta} = 1 - 0.0033 \\mid \\beta \\mid

    with a maximum of :math:`\\gamma_{\\beta} = 0.736` for
    :math:`\\mid \\beta \\mid > 80`

    Parameters
    ----------
    beta : float
        the angle of wave attack [rad]

    Returns
    -------
    gamma_beta : float
        the influence factor for oblique wave attack
    """
    beta_deg = beta*180/np.pi

    if beta_deg <= 80:
        return 1 - 0.0033 * np.abs(beta_deg)
    else:
        return 0.736

def rubble_mound(
        Hm0, q, xi_m_min_1, alpha, beta, gamma_b, gamma_v, gam_star,
        armour_layer, layers=1, permeability='permeable', safety=1,
        Gc=None, Dn50=None, limit=True):
    """ Compute the crest freeboard of a rubble mound breakwater

    Computes the crest freeboard of a rubble mound breakwater using
    equation 5.10 and 5.12 from EurOtop (2018).

    .. math::
       \\frac{q}{\\sqrt{g \\cdot H_{m 0}^{3}}}=\\frac{0.023}{\\sqrt{\\tan
       \\alpha}}\\gamma_{b} \\cdot \\xi_{m-1,0} \\cdot \\exp \\left[-
       \\left(2.7 \\frac{R_{c}}{\\xi_{m-1,0} \\cdot H_{m 0} \\cdot
       \\gamma_{b} \\cdot \\gamma_{f} \\cdot \\gamma_{\\beta} \\cdot
       \\gamma_{v}}\\right)^{1.3}\\right]

    with a maximum of

    .. math::
       \\frac{q}{\\sqrt{g \\cdot H_{m 0}^{3}}}=0.09 \\cdot \\exp \\left
       [-\\left(1.5 \\frac{R_{c}}{H_{m 0} \\cdot \\gamma_{f} \\cdot
       \\gamma_{\\beta} \\cdot \\gamma^{*}}\\right)^{1.3}\\right]

    The reliability of the first equation is given by
    :math:`\\sigma` (0.023) = 0.003 and :math:`\\sigma` (2.7) = 0.20,
    and for the second equation by :math:`\\sigma` (0.09) = 0.0135 and
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
    beta : float
        the angle of wave attack [rad]
    gamma_b : float
        influence factor for a berm [-]
    gamma_v : float
        influence factor for a vertical wall [-]
    gamma_star : float
        overall influence factor for a storm wall on slope or promenade [-]
    armour_layer : str
        name of the material of the armour layer, supported materials:
        Smooth, Rock, Cubes, Antifers, HARO, Tetrapods, Dolose,
        Accropode I, Xbloc, CoreLoc, Accropode II, Cubipods
    layers : {1, 2}, default: 1
        number of layers in the armour layer, required if the armour
        layer is made out of: Rock, Cubes or Cubipods
    permeability : {'permeable', 'impermeable'}, default: 'permeable'
        permeability of the armour layer, required if the armour layer
        is made out of Rock
    safety : float, optional, default: 1
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
    g = 9.81

    gam_f = gamma_f(
        armour_layer=armour_layer, xi_m_min_1=xi_m_min_1, layers=layers,
        permeability=permeability)
    gam_beta = gamma_beta(beta)

    constant_a = 0.023 + safety*0.003
    constant_b = 2.7 - safety*0.2
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
    arg = (q * np.sqrt(np.tan(alpha))
           / (constant_a*gamma_b*xi_m_min_1*np.sqrt(g*Hm0**3)))
    arg_upper = q/(constant_c*np.sqrt(g*Hm0**3))

    if arg <= 1:
        Rc = ((xi_m_min_1*Hm0*gamma_b*gam_f*gam_beta*gamma_v/constant_b)
              * (-np.log(arg))**(1/1.3))
    else:
        Rc = 0
        user_warning('Encountered negative freeboard')

    if arg_upper <= 1: # handle for RuntimeWarning
        Rc_upper = ((Hm0*gam_f*gam_beta*gam_star/constant_d)
                    * (-np.log(arg_upper))**(1/1.3))
    else:
        Rc_upper = 0
        user_warning('Encountered negative freeboard in upper bound')

    Rc = min(Rc, Rc_upper)

    return Rc

def vertical_deep(Hm0, q, safety=1, limit=True):
    """ Rc if the foreshore does not have an influence

    Compute the crest freeboard for a vertical or composite vertical
    wall if the foreshore does not have an influence. Implementation
    of equation 7.1 from EurOtop (2018).

    .. math::
       \\frac{q}{\\sqrt{g \\cdot H_{m 0}^{3}}}=0.047 \\cdot \\exp \\left[
       -\\left(2.35 \\frac{R_{c}}{H_{m 0}}\\right)^{1.3}\\right]

    The reliability of the equation is given by
    :math:`\\sigma` (0.047) = 0.007 and :math:`\\sigma` (2.35) = 0.23.

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    safety : float, optional, default: 1
        safety factor of the design, positive values increase the safety
        of the design by increasing the mean value of the model constants
        with the number of standard deviations specified. In accordance
        with the recommendation from EurOtop (2018), the default value is
        set to 1 standard deviation.
    limit : bool, optional, default: True
        If True, the discharge will be set to the limit for zero
        discharge in case the given discharge is below this limit. If
        False, the discharge will not be changed.

    Returns
    -------
    Rc : float
        the crest freeboard of the structure [m]
    """
    g = 9.81

    # check given overtop discharge against EurOtop limit for zero discharge
    limit_zero_q = np.sqrt(g*Hm0**3) * 0.01
    if limit_zero_q > q:
        msg = (f'Given value of q = {q}, is below EurOtop limit for zero '
               f'discharge, which is {np.round(limit_zero_q, 3)}')

        if limit:
            msg += f'. Using {np.round(limit_zero_q, 3)} instead as limit=True'
            q = limit_zero_q

        user_warning(msg)

    # convert q from l to m^3
    q = q/1000

    constant_a = 0.047 + safety*0.007
    constant_b = 2.35 - safety*0.23

    Rc = (1/constant_b*Hm0
          * (-np.log(q/(constant_a*np.sqrt(g*Hm0**3))))**(1/1.3))
    return Rc

def vertical_no_breaking(Hm0, q, safety=1, limit=True):
    """ Rc if no possibility for breaking waves

    Compute the crest freeboard for a vertical or composite vertical
    wall if there are no breaking waves. Implementation of equation 7.5
    from EurOtop (2018)

    .. math::
       \\frac{q}{\\sqrt{g H_{m 0}^{3}}}=0.05 \\exp \\left(-2.78
       \\frac{R_{c}}{H_{m 0}}\\right)

    The reliability of the equation is given by
    :math:`\\sigma` (0.05) = 0.012 and :math:`\\sigma` (2.78) = 0.17.

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    safety : float, optional, default: 1
        safety factor of the design, positive values increase the safety
        of the design by increasing the mean value of the model constants
        with the number of standard deviations specified. In accordance
        with the recommendation from EurOtop (2018), the default value is
        set to 1 standard deviation.
    limit : bool, optional, default: True
        If True, the discharge will be set to the limit for zero
        discharge in case the given discharge is below this limit. If
        False, the discharge will not be changed.

    Returns
    -------
    Rc : float
        the crest freeboard of the structure [m]
    """
    g = 9.81

    # check given overtop discharge against EurOtop limit for zero discharge
    limit_zero_q = np.sqrt(g*Hm0**3) * 0.01
    if limit_zero_q > q:
        msg = (f'Given value of q = {q}, is below EurOtop limit for zero '
               f'discharge, which is {np.round(limit_zero_q, 3)}')

        if limit:
            msg += f'. Using {np.round(limit_zero_q, 3)} instead as limit=True'
            q = limit_zero_q

        user_warning(msg)

    # convert q from l to m^3
    q = q/1000

    constant_a = 0.05 + safety*0.012
    constant_b = 2.78 - safety*0.17

    Rc = 1/constant_b*Hm0 * -np.log(q/(constant_a*np.sqrt(g*Hm0**3)))
    return Rc

def vertical_normal(Hm0, q, h, s_m_min_1, safety=1, limit=True):
    """ Rc for a vertical wall if normal freeboard and breaking waves

    Compute the crest freeboard for a vertical wall if there is a
    possibility for breaking waves, but the freeboard is not low.
    Implementation of equation 7.8 from EurOtop (2018)

    .. math::
       \\frac{q}{\\sqrt{g H_{m 0}^{3}}}=0.0014\\left(\\frac{H_{m 0}}
       {h s_{m-1,0}}\\right)^{0.5}\\left(\\frac{R_{c}}{H_{m 0}}
       \\right)^{-3}

    The reliability of the equation is given by
    :math:`\\sigma` (0.0014) = 0.0006

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    h : float
        water depth in front of the toe of the structure [m]
    s_m_min_1 : float
        :math:`s_{m-1.0}`, wave steepness with the spectral wave length
        (:math:`L_{m-1.0}`) [-]
    safety : float, optional, default: 1
        safety factor of the design, positive values increase the safety
        of the design by increasing the mean value of the model constants
        with the number of standard deviations specified. In accordance
        with the recommendation from EurOtop (2018), the default value is
        set to 1 standard deviation.
    limit : bool, optional, default: True
        If True, the discharge will be set to the limit for zero
        discharge in case the given discharge is below this limit. If
        False, the discharge will not be changed.

    Returns
    -------
    Rc : float
        the crest freeboard of the structure [m]
    """
    g = 9.81

    # check given overtop discharge against EurOtop limit for zero discharge
    limit_zero_q = np.sqrt(g*Hm0**3) * 0.01
    if limit_zero_q > q:
        msg = (f'Given value of q = {q}, is below EurOtop limit for zero '
               f'discharge, which is {np.round(limit_zero_q, 3)}')

        if limit:
            msg += f'. Using {np.round(limit_zero_q, 3)} instead as limit=True'
            q = limit_zero_q

        user_warning(msg)

    # convert q from l to m^3
    q = q/1000

    constant = 0.0014 + safety*0.0006

    Rc = Hm0 * (q/(np.sqrt(Hm0*g*Hm0**3/(h*s_m_min_1))*constant))**(1/-3)
    return Rc

def vertical_low(Hm0, q, h, s_m_min_1, safety=1, limit=True):
    """ Rc for a vertical wall if low freeboard and breaking waves

    Compute the crest freeboard for a vertical wall if there is a
    possibility for breaking waves, and the freeboard is low.
    Implementation of equation 7.7 from EurOtop (2018)

    .. math::
       \\frac{q}{\\sqrt{g H_{m 0}^{3}}}=0.011\\left(\\frac{H_{m 0}}{h
       s_{m-1,0}}\\right)^{0.5} \\exp \\left(-2.2 \\frac{R_{c}}{H_{m 0}}
       \\right)

    The reliability of the equation is given by
    :math:`\\sigma` (0.011) = 0.0045

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    h : float
        water depth in front of the toe of the structure [m]
    s_m_min_1 : float
        :math:`s_{m-1.0}`, wave steepness with the spectral wave length
        (:math:`L_{m-1.0}`) [-]
    safety : float, optional, default: 1
        safety factor of the design, positive values increase the safety
        of the design by increasing the mean value of the model constants
        with the number of standard deviations specified. In accordance
        with the recommendation from EurOtop (2018), the default value is
        set to 1 standard deviation.
    limit : bool, optional, default: True
        If True, the discharge will be set to the limit for zero
        discharge in case the given discharge is below this limit. If
        False, the discharge will not be changed.

    Returns
    -------
    Rc : float
        the crest freeboard of the structure [m]
    """
    g = 9.81

    # check given overtop discharge against EurOtop limit for zero discharge
    limit_zero_q = np.sqrt(g*Hm0**3) * 0.01
    if limit_zero_q > q:
        msg = (f'Given value of q = {q}, is below EurOtop limit for zero '
               f'discharge, which is {np.round(limit_zero_q, 3)}')

        if limit:
            msg += f'. Using {np.round(limit_zero_q, 3)} instead as limit=True'
            q = limit_zero_q

        user_warning(msg)

    # convert q from l to m^3
    q = q/1000

    constant = 0.011 + safety*0.0045

    Rc = Hm0/-2.2 * np.log(q/(np.sqrt(Hm0*g*Hm0**3/(h*s_m_min_1))*constant))
    return Rc

def composite_normal(Hm0, q, h, d, s_m_min_1, safety=1, limit=True):
    """ Rc for a composite wall if normal freeboard and breaking waves

    Compute the crest freeboard for a composite vertical wall if there
    is a possibility for breaking waves, but the freeboard is not low.
    Implementation of equation 7.14 from EurOtop (2018)

    .. math::
       \\frac{q}{\\sqrt{g H_{m 0}^{3}}}=1.3\\left(\\frac{d}{h}\\right)^{0.5}
       0.0014\\left(\\frac{H_{m 0}}{h s_{m-1,0}}\\right)^{0.5}\\left(
       \\frac{R_{c}}{H_{m 0}}\\right)^{-3}

    The reliability of the equation is given by
    :math:`\\sigma` (0.0014) = 0.0006.

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    h : float
        water depth in front of the toe of the structure [m]
    d : float
        water depth above the toe of the structure [m]
    s_m_min_1 : float
        :math:`s_{m-1.0}`, wave steepness with the spectral wave length
        (:math:`L_{m-1.0}`) [-]
    safety : float, optional, default: 1
        safety factor of the design, positive values increase the safety
        of the design by increasing the mean value of the model constants
        with the number of standard deviations specified. In accordance
        with the recommendation from EurOtop (2018), the default value is
        set to 1 standard deviation.
    limit : bool, optional, default: True
        If True, the discharge will be set to the limit for zero
        discharge in case the given discharge is below this limit. If
        False, the discharge will not be changed.

    Returns
    -------
    Rc : float
        the crest freeboard of the structure [m]
    """
    g = 9.81

    # check given overtop discharge against EurOtop limit for zero discharge
    limit_zero_q = np.sqrt(g*Hm0**3) * 0.01
    if limit_zero_q > q:
        msg = (f'Given value of q = {q}, is below EurOtop limit for zero '
               f'discharge, which is {np.round(limit_zero_q, 3)}')

        if limit:
            msg += f'. Using {np.round(limit_zero_q, 3)} instead as limit=True'
            q = limit_zero_q

        user_warning(msg)

    # convert q from l to m^3
    q = q/1000

    constant = 0.0014 + safety*0.0006

    Rc = Hm0 * (q/(np.sqrt(d*g*Hm0**4/(h**2*s_m_min_1))*1.3*constant))**(1/-3)
    return Rc

def composite_low(Hm0, q, h, d, s_m_min_1, safety=1, limit=True):
    """ Rc for a composite wall if low freeboard and breaking waves

    Compute the crest freeboard for a composite vertical wall if there
    is a possibility for breaking waves, and the freeboard is low.
    Implementation of equation 7.15 from EurOtop (2018)

    .. math::
       \\frac{q}{\\sqrt{g H_{m 0}^{3}}}=1.3\\left(\\frac{d}{h}
       \\right)^{0.5} 0.011\\left(\\frac{H_{m 0}}{h s_{m-1,0}}
       \\right)^{0.5} \\exp \\left(-2.2 \\frac{R_{c}}{H_{m 0}}\\right)

    The reliability of the equation is given by
    :math:`\\sigma` (0.011) = 0.0045.

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    h : float
        water depth in front of the toe of the structure [m]
    d : float
        water depth above the toe of the structure [m]
    s_m_min_1 : float
        :math:`s_{m-1.0}`, wave steepness with the spectral wave length
        (:math:`L_{m-1.0}`) [-]
    safety : float, optional, default: 1
        safety factor of the design, positive values increase the safety
        of the design by increasing the mean value of the model constants
        with the number of standard deviations specified. In accordance
        with the recommendation from EurOtop (2018), the default value is
        set to 1 standard deviation.
    limit : bool, optional, default: True
        If True, the discharge will be set to the limit for zero
        discharge in case the given discharge is below this limit. If
        False, the discharge will not be changed.

    Returns
    -------
    Rc : float
        the crest freeboard of the structure [m]
    """
    g = 9.81

    # check given overtop discharge against EurOtop limit for zero discharge
    limit_zero_q = np.sqrt(g*Hm0**3) * 0.01
    if limit_zero_q > q:
        msg = (f'Given value of q = {q}, is below EurOtop limit for zero '
               f'discharge, which is {np.round(limit_zero_q, 3)}')

        if limit:
            msg += f'. Using {np.round(limit_zero_q, 3)} instead as limit=True'
            q = limit_zero_q

        user_warning(msg)

    # convert q from l to m^3
    q = q/1000

    constant = 0.011 + safety*0.0045
    Rc = Hm0/-2.2*np.log(q/(np.sqrt(d*g*Hm0**4/(h**2*s_m_min_1))*1.3*constant))
    return Rc

def vertical(
        Hm0, q, h, d, L_m_min_1, s_m_min_1, safety=1, logger=None, limit=True):
    """ Compute crest freeboard for vertical and composite vertical walls

    Compute the crest freeboard, Rc, of a vertical or composite vertical
    wall for a given mean overtop discharge. The function is an
    implementation of the decision chart, figure 7.2, from EurOtop
    (2018). The function determines if the input classifies as a
    vertical or composite vertical wall, if breaking is possible and if
    the structure has a low freeboard. Based on the classification the
    function calls the corresponding function and computes the freeboard.

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    h : float
        water depth in front of the toe of the structure [m]
    d : float
        water depth above the toe of the structure [m]
    L_m_min_1 : float
        :math:`L_{m-1.0}`, spectral wave length in deep water [m]
    s_m_min_1 : float
        :math:`s_{m-1.0}`, wave steepness with the spectral wave length
        (:math:`L_{m-1.0}`) [-]
    safety : float, optional, default: 1
        safety factor of the design, positive values increase the safety
        of the design by increasing the mean value of the model constants
        with the number of standard deviations specified. In accordance
        with the recommendation from EurOtop (2018), the default value is
        set to 1 standard deviation.
    logger : dict, optional, default: None
        dict to log messages, must have keys 'INFO' and 'WARNINGS'
    limit : bool, optional, default: True
        If True, the discharge will be set to the limit for zero
        discharge in case the given discharge is below this limit. If
        False, the discharge will not be changed.

    Returns
    -------
    Rc : float
        the crest freeboard of the structure [m]
    """
    if h/Hm0 > 4:
        Rc = vertical_deep(Hm0, q, safety, limit)
        log = 'overtopping for vertical in deep water'
    else:
        if d/h > 0.6:  # check if vertical wall
            if h**2/(Hm0*L_m_min_1) < 0.23:  # check if waves are breaking
                Rc = vertical_normal(Hm0, q, h, s_m_min_1, safety, limit)
                log = 'overtopping for vertical with normal freeboard'
                if Rc/Hm0 < 1.35:  # check if freeboard is low
                    Rc = vertical_low(Hm0, q, h, s_m_min_1, safety, limit)
                    log = 'overtopping for vertical with low freeboard'
            else:
                Rc = vertical_no_breaking(Hm0, q, safety, limit)
                log = 'overtopping for vertical with no breaking waves'
        else:  # not a vertical wall, treat as composite wall
            if h*d/(Hm0*L_m_min_1) < 0.65:  # check if waves are breaking
                Rc = composite_normal(Hm0, q, h, d, s_m_min_1, safety, limit)
                log = 'overtopping for composite with normal freeboard'
                if Rc/Hm0 < 1.35:  # check if freeboard is low
                    Rc = composite_low(Hm0, q, h, d, s_m_min_1, safety, limit)
                    log = 'overtopping for composite with low freeboard'
            else:
                Rc = vertical_no_breaking(Hm0, q, safety, limit)
                log = 'overtopping for composite with no breaking waves'
    if logger != None:
        logger['INFO'].append(log)

    return Rc
