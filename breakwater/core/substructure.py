from ..utils.exceptions import user_warning, InputError, NotSupportedError

# supported armour layers, i.e. for these armour layers the rules for
# the underlayer and filter have been implemented
def _supported_armour_layers():
    """ returns the supported types of material for the armour layer """
    return ['Rock', 'Xbloc', 'XblocPlus']


def underlayer(Dn_armour, armour_layer, rho, rho_rock=2650):
    """ Design first underlayer of a rubble mound breakwater

    Design the layer of rock immediately below the armour layer, the
    nominal diameter of this layer (the underlayer) is determined using
    one of the following rules:

    - | Rock (CIRIA, CUR, CETMEF, 2007, p. 630)
      | :math:`M_{50u} = \\frac{M_{50a}}{15}` to :math:`\\frac{M_{50a}}{10}`

    - | Xbloc (Delta Marine Consultants, 2018)
      | :math:`M_{50u} = \\frac{M_{50a}}{15}` to :math:`\\frac{M_{50a}}{6}`

    - | XblocPlus (Delta Marine Consultants, 2018)
      | :math:`M_{50u} = \\frac{M_{50a}}{20}` to :math:`\\frac{M_{50a}}{8}`

    Parameters
    ----------
    Dn_armour : float
        nominal diameter of the armour layer [m]
    armour_layer : {'Rock', 'Xbloc', 'XblocPlus'}
        type of the armour layer
    rho : float, optional, default: 2650
        density of the armourstone [kg/m³]
    rho_rock : float, optional, default: 2650
        density of the rock used in the underlayer [kg/m³]

    Returns
    -------
    [float, float]
        list with the upper bound and lower bound limit of the Dn50 of
        the underlayer
    """
    M_armour = rho*Dn_armour**3

    if armour_layer == 'Rock':
        M_upper = M_armour/10
        M_lower = M_armour/15

    elif armour_layer == 'Xbloc':
        M_upper = M_armour/6
        M_lower = M_armour/15

    elif armour_layer == 'XblocPlus':
        M_upper = M_armour/8
        M_lower = M_armour/20

    else:
        supported = ', '.join(_supported_armour_layers())
        raise NotSupportedError(
            (f'{armour_layer} is not implemented as a rule for the underlayer'
             f', supported types are: {supported}.'))

    dn_upper = (M_upper/rho_rock)**(1/3)
    dn_lower = (M_lower/rho_rock)**(1/3)

    return [dn_lower, dn_upper]

def filter_layers(Dn, rho=2650):
    """ Design filter layer of a rubble mound breakwater

    Design a layer of rock below the first underlayer, or subsequent
    layers of rock. The nominal diameter of this layer (the filter layer)
    is given by the following rule:

    - | Rock (Van den Bos and Verhagen, 2018, p. 142)
      | :math:`M_{50l} = \\frac{M_{50u}}{25}` to :math:`\\frac{M_{50u}}{10}`

    Parameters
    ----------
    Dn : float
        nominal diameter of the layer above the filter layer [m]
    rho : float, optional, default: 2650
        density of the armourstone [kg/m³]

    Returns
    -------
    [float, float]
        list with the upper bound and lower bound limit of the Dn50 of
        the filter layer
    """
    M_upper_layer = Dn**3

    dn_upper = Dn/(10**(1/3))
    dn_lower = Dn/(25**(1/3))

    return [dn_lower, dn_upper]

def layer_coefficient(material, layers=None, placement=None):
    """ Get the layer thickness coefficient

    Determine the layer thickness coefficient, :math:`k_{t}`. In Table
    1 and 2 the included layer thickness coefficients can be found.

    Table 1: :math:`k_{t}` values for Rock (CIRIA, CUR, CETMEF, 2007, p. 126)

    +---------+-----------+-------+
    | Layers  | Placement | Value |
    +=========+===========+=======+
    | 1       | Dense     | 0.84  |
    +---------+-----------+-------+
    | 2       | Standard  | 0.91  |
    +---------+-----------+-------+
    | 2       | Dense     | 0.91  |
    +---------+-----------+-------+

    Table 2: :math:`k_{t}` values for Armour Units (CIRIA, CUR, CETMEF, 2007, p. 260)

    +-------------+--------+-------+
    | Armour Unit | Layers | Value |
    +=============+========+=======+
    | Cubes       |   2    |  1.1  |
    +-------------+--------+-------+
    | Tetrapods   |   2    |  1.02 |
    +-------------+--------+-------+
    | Dolos       |   2    |  0.94 |
    +-------------+--------+-------+
    | Accropode   |   1    |  1.29 |
    +-------------+--------+-------+
    | CoreLoc     |   1    |  1.52 |
    +-------------+--------+-------+
    | Xbloc       |   1    |  1.4  |
    +-------------+--------+-------+
    | XblocPlus   |   1    |  1.33 |
    +-------------+--------+-------+

    Parameters
    ----------
    material : str
        material of the layer
    layers : int, optional, default: None
        number of layers, required if the material is rock. For armour
        units the number of layers is not required, but a warning will
        be shown if the specified layer is different from the table
    placement : {Standard, Dense}, optional, default: None
        placement of the material, required if the material is rock

    Returns
    -------
    kt : float
        the layer thickness coefficient
    """
    # define dict for rock
    rock = {

        1: {'dense': 0.84},

        2: {'standard': 0.91, 'dense': 0.91}

    }

    # define dict for armour units
    units = {

        'Cubes': {'kt': 1.1, 'layers': 2},

        'Tetrapods': {'kt': 1.02, 'layers': 2},

        'Dolos': {'kt': 0.94, 'layers': 2},

        'Accropode': {'kt': 1.29, 'layers': 1},

        'CoreLoc': {'kt': 1.52, 'layers': 1},

        'Xbloc': {'kt': 1.4, 'layers': 1},

        'XblocPlus': {'kt': 1.33, 'layers': 1}

    }
    # Check if material has Accropode in the name,
    # as it can be Accropode I or Accropode II
    if 'Accropode' in material:
        # set material to Accropode
        material = 'Accropode'

    # get the kt value
    if material is 'Rock':
        # check if layer is given
        if layers is None:
            raise InputError(
                'Layers is a required argument if the material is Rock')

        # check if placement is given
        if placement is None:
            raise InputError(
                'Placement is a required argument if the material is Rock')

        # get the kt value
        kt = rock[layers][placement.lower()]

    elif material in list(units.keys()):
        # get kt value
        kt = units[material]['kt']

        # check layers
        num_layers = units[material]['layers']
        if layers is not None and layers != num_layers:
            # check number of layers
            if num_layers == 1:
                msg = (f'kt value for {material} is only determined for '
                       f'{num_layers} layer, not {layers}')
            else:
                # plural
                msg = (f'kt value for {material} is only determined for '
                       f'{num_layers} layers, not {layers}')

            user_warning(msg)

    else:
        # not supported armour layer
        kt = 1
        user_warning(f'{material} is not supported, continued with kt=1')

    return kt
