import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pickle
from warnings import catch_warnings
from tabulate import tabulate

from .utils._kwarg_validator import _process_kwargs, _RM_vkwargs, _C_vkwargs
from .utils._design_explorer import _DE_params
from .utils.exceptions import RockGradingError, ArmourUnitsError, InputError, user_warning, NotSupportedError
from .utils._progress import ProgressBar
from .utils._excel_utils import _convert_string
from .utils.cost import _process_cost, cost_influence
from .caisson import Caisson
from .conditions import LimitState
from .material import read_grading, read_units
from .rubble import RockRubbleMound, ConcreteRubbleMound, ConcreteRubbleMoundRevetment


def read_configurations(
        filepath, structure, kd=None, name=None, LS=None, Grading=None,
        ArmourUnits=None, BermMaterial=None):
    """ Conceptual design for multiple breakwaters from an Excel file

    Make a conceptual design for multiple (types) of breakwaters from
    an Excel input file. The Excel input file can be generated with
    :obj:`bw.generate_excel <breakwater.utils.input_generator.generate_excel>`

    .. warning:
       When reading Xbloc or XblocPlus armour units from an Excel file,
       the correction factor applied to the nominal diameter is not
       computed. Specify :py:class:`Xbloc` or :py:class:`XblocPlus`
       for ArmourUnits so that the most optimal design is made.

    Parameters
    ----------
    structure : {'RRM', 'CRM', 'CRMR', 'RC', 'CC'}
        structure for which conceptual designs must be generated. RRM
        for a rubble mound with rock as armour layer, CRM for a rubble
        mound with concrete armour units as armour layer, CRMR for a rubble
        mound revetment with concrete armour units as armour layer, RC for a
        vertical (composite) breakwater with rock as armour layer for
        the foundation and CC for a vertical (composite) breakwater with
        concrete armour units as armour layer for the foundation.
    kd : int, optional, default: None
        stability coefficient [-]
    name : str, optional, default: None
        name of the ArmourUnit
    LS : py:class:`LimitState`, optional, default: None
        ULS, SLS or another limit state defined with
        :py:class:`LimitState`, by default the LimitState is read from
        the Excel input file
    Grading : :py:class:`RockGrading`, optional, default: None
        standard rock grading defined in the NEN-EN 13383-1 or a user
        defined rock grading. By default the Grading is read from
        the Excel input file
    ArmourUnit : obj, optional, default: None
        armour unit class which inherits from :py:class:`ConcreteArmour`,
        for instance :py:class:`Xbloc` or :py:class:`XblocPlus`. By
        default the ArmourUnit is read from the Excel input file. This
        argument is used for RRM.
    BermMaterial : obj, optional, default: None
        armour unit class which inherits from :py:class:`ConcreteArmour`,
        for instance :py:class:`Xbloc` or :py:class:`XblocPlus`. By
        default the BermMaterial is read from the Excel input file.
        This argument is used for CC.

    Returns
    -------
    bw.Configurations
        a Configurations object with breakwater concepts
    """
    # convert the input of structure to a list
    if isinstance(structure, list):
        # must be a list so no change
        structure = structure
    elif isinstance(structure, str):
        # convert single input to list
        structure = [structure]

    # import the excel as a DataFrame
    df = pd.read_excel(filepath, skiprows=1, sheet_name='Parameters')

    # generate input dict
    params = {}

    # load LimitState
    if LS is None:
        # load LimitState data from excel file
        df_ls = pd.read_excel(filepath, index_col=0, sheet_name='LimitState')

        # convert to dict
        input_ls = df_ls.to_dict()['Value']

        # generate LimitState object
        params['LimitState'] = LimitState(**input_ls)
    else:
        params['LimitState'] = LS

    # load RockGrading
    if Grading is None:
        # load Grading
        params['Grading'] = read_grading(filepath, sheet_name='RockGrading')
    else:
        params['Grading'] = Grading

    # load ArmourUnits
    if ArmourUnits is None and 'CRM' in structure:
        # check if kd and name are given
        if kd is None:
            raise InputError(
                ('Missing argument: kd, when reading ArmourUnits from Excel '
                 'the kd value must be given as input'))
        if name is None:
            raise InputError(
                ('Missing argument: name, when reading ArmourUnits from Excel '
                 'the name of the ArmourUnit must be given as input'))

        # load ArmourUnits
        params['ArmourUnit'] = read_units(
            filepath, kd=kd, name=name, sheet_name='ArmourUnit')
    elif ArmourUnits is not None and 'CRM' in structure:
        params['ArmourUnit'] = ArmourUnit

    # load BermMaterial
    if BermMaterial is None and 'CC' in structure:
        # check if kd and name are given
        if kd is None:
            raise InputError(
                ('Missing argument: kd, when reading ArmourUnits from Excel '
                 'the kd value must be given as input'))
        if name is None:
            raise InputError(
                ('Missing argument: name, when reading ArmourUnits from Excel '
                 'the name of the ArmourUnit must be given as input'))

        # load ArmourUnits
        params['BermMaterial'] = read_units(
            filepath, kd=kd, name=name, sheet_name='BermMaterial')
    elif BermMaterial is not None and 'CC' in structure:
        params['BermMaterial'] = BermMaterial

    # set column names and drop validation
    df.columns = ['parameter', 'value', 'min', 'max', 'num', 'validation']
    df.drop(labels=['validation'],axis=1, inplace=True)

    parameters = df['parameter'].values
    values = df['value'].values

    for i, parameter in enumerate(parameters):
        if 'slope' in parameter:
            # check if value is nan
            if not isinstance(values[i], str) and np.isnan(values[i]):
                # parameter is varying
                # convert min and max to tuple
                min = _convert_string(df['min'].values[i], 'slope', parameter)
                max = _convert_string(df['max'].values[i], 'slope', parameter)

                # get num
                num = df['num'].values[i]

                # set value
                value = (min, max, int(num))

            else:
                # parameter is set as constant parameter
                value = _convert_string(values[i], 'slope', parameter)

        elif 'lambda' in parameter:
            # convert lambda to list
            value = _convert_string(values[i], 'lambda', parameter)

        else:
            # get the value
            value = values[i]

            # check if filter rule
            if isinstance(value, str):
                params[parameter] = value
                continue

            # check if value is nan
            if np.isnan(value):
                # get varying input
                min = df['min'].values[i]
                max = df['max'].values[i]
                num = df['num'].values[i]

                # check if one of the values is nan
                if np.isnan([min, max, num]).all():
                    # if all are nan set value to None
                    value = None
                    continue

                elif np.isnan([min, max, num]).any():
                    raise TypeError(
                        (f'Varying parameter {parameter} is incorrectly '
                          'formatted'))

                # convert to tuple
                value = (min, max, int(num))

        # add to params dict
        params[parameter] = value

    return Configurations(structure=structure, **params)

def read_breakwaters(filepath):
    """ Read a breakwaters file into Configurations

    Parameters
    ----------
    filepath : path
        any valid string path is acceptable, filepath must end with
        .breakwaters

    Returns
    -------
    bw.Configurations
        a Configurations object with breakwater concepts
    """
    # get the extension of the file
    extension = os.path.splitext(filepath)[1]

    # check if a extension has been included in the filepath
    if not extension:
        # if not add .breakwaters to the filepath
        filepath = f'{filepath}.breakwaters '
    else:
        # check if extension is .breakwaters
        if extension == '.breakwaters':
            pass
        else:
            raise InputError(
                (f'Can\'t load Configurations from {extension}, must be '
                  '.breakwaters'))

    # load breakwaters file, and return Configurations object
    file = open(filepath, 'rb')
    config = pickle.load(file)

    return config


class Configurations:
    """ Make a conceptual design for multiple breakwaters

    With this class a parametric design of one or more types of
    breakwaters. The types of structures included in the tool are: a
    rubble mound breakwater with rock (RRM), concrete armour units
    (CRM) as armour layer, concrete armour units as armour layer as 
    part of a revetment  and a vertical (composite) breakwater with 
    rock (RC) or concrete armour units (CC) as armour layer of the
    foundation layer. Each structure has its own set of required
    (table 1) and optional (table 2) parameters.

    For the parametric design some parameters are allowed to vary,
    these parameters have a superscript in table 1 and 2. The input of
    a varying parameter must be a tuple with a length of 3, i.e.
    (min, max, num), where min and max are the minimum and maximum
    value and num is the number of samples to generate.

    .. note::
       Alternatively, a design can also be made from an Excel input file
       by using :py:func:`read_configurations`. The required Excel input
       file can be generated by using
       :obj:`bw.generate_excel <breakwater.utils.input_generator.generate_excel>`

    Table 1: required parameters

    +--------------------+-----+-----+-----+-----+-----+
    | Parameter          | RRM | CRM | CRMR| RC  | CC  |
    +====================+=====+=====+=====+=====+=====+
    | LimitState         |  x  |  x  |  x  |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | rho_w              |  x  |  x  |  x  |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | slope_foreshore    |  x  |  x  |  x  |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | Grading            |  x  |  x  |  x  |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | slope :sup:`1`     |  x  |  x  |  x  |     |     |
    +--------------------+-----+-----+-----+-----+-----+
    | B :sup:`1`         |  x  |  x  |  x  |     |     |
    +--------------------+-----+-----+-----+-----+-----+
    | Dn50_core :sup:`1` |  x  |  x  |  x  |     |     |
    +--------------------+-----+-----+-----+-----+-----+
    | N                  |  x  |     |     |     |     |
    +--------------------+-----+-----+-----+-----+-----+
    | ArmourUnit         |     |  x  |  x  |     |     |
    +--------------------+-----+-----+-----+-----+-----+
    | Pc :sup:`1`        |     |     |     |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | rho_c              |     |     |     |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | rho_fill           |     |     |     |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | Bm :sup:`1`        |     |     |     |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | hb :sup:`1`        |     |     |     |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | mu                 |     |     |     |  x  |  x  |
    +--------------------+-----+-----+-----+-----+-----+
    | BermMaterial       |     |     |     |     |  x  |
    +--------------------+-----+-----+-----+-----+-----+

    :sup:`1` Parameter is allowed to vary, enter as a tuple with
    (min, max, num), where min and max are the minimum and maximum
    value and num is the number of samples to generate.

    Table 2: Optional parameters with default values

    +---------------------------+---------+-----+-----+-----+-----+-----+
    | Parameter                 | Default | RRM | CRM | CRMR| RC  | CC  |
    +===========================+=========+=====+=====+=====+=====+=====+
    | safety                    |    1    |  o  |  o  |  o  |  o  |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | beta                      |    0    |  o  |  o  |  o  |  o  |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | Soil                      |  None   |  o  |  o  |  o  |  o  |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | slope_toe :sup:`1`        |  (2,3)  |  o  |  o  |  o  |     |     |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | phi                       |    40   |  o  |  o  |  o  |     |     |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | B_toe :sup:`1`            |  None   |  o  |  o  |  o  |     |     |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | vdm                       |   max   |  o  |     |     |     |     |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | layers_rock               |    2    |  o  |     |     |  o  |     |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | layers_units              |    1    |     |  o  |  o  |     |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | layers_underlayer         |    2    |  o  |  o  |  o  |     |     |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | filter_rule               |  None   |     |  o  |  o  |     |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | pe_max                    |   500   |     |     |     |  o  |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | SF_sliding                |   1.2   |     |     |     |  o  |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | SF_turning                |   1.2   |     |     |     |  o  |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | slope_foundation :sup:`1` |  (2,3)  |     |     |     |  o  |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+
    | \lambda_                  | [1,1,1] |     |     |     |  o  |  o  |
    +---------------------------+---------+-----+-----+-----+-----+-----+


    :sup:`1` Parameter is allowed to vary, enter as a tuple with
    (min, max, num), where min and max are the minimum and maximum
    value and num is the number of samples to generate.

    Parameters
    ----------
    structure : {'RRM', 'CRM', 'CRMR', 'RC', 'CC'}
        structure for which conceptual designs must be generated. RRM
        for a rubble mound with rock as armour layer, CRM for a rubble
        mound with concrete armour units as armour layer, CRMR for a rubble
        mound revetment with concrete armour units as armour layer. RC for a
        vertical (composite) breakwater with rock as armour layer for
        the foundation and CC for a vertical (composite) breakwater with
        concrete armour units as armour layer for the foundation.
    LimitState : :py:class:`LimitState` or list of :py:class:`LimitState`
        ULS, SLS or another limit state defined with
        :py:class:`LimitState`
    rho_w : float
        density of water [kg/m³]
    slope_foreshore : tuple
        slope of the foreshore (V, H). For example a slope of 1:100 is
        defined as (1, 100)
    Grading : :py:class:`RockGrading`
        standard rock grading defined in the NEN-EN 13383-1 or a user
        defined rock grading.
    safety : float, optional, default: 1
        safety factor of design (number of standard deviations from the
        mean)

    Keyword arguments
    -----------------
    slope : tuple
        Slope of the armour layer (V,H). For example a slope of 3V:4H
        is defined as (3,4). Required for RRM and CRM.
    B : float
        Crest width [m], required for RRM and CRM
    Dn50_core : float
        nominal diameter for the stones in the core of the breakwater [m],
        required for RRM and CRM
    N : int
        Number of incident waves at the toe of the structure [-],
        required for RRM
    ArmourUnit : obj
        armour unit class which inherits from :py:class:`ConcreteArmour`,
        for instance :py:class:`Xbloc` or :py:class:`XblocPlus`. Required
        argument for CRM.
    Pc : float
        contribution of concrete to the total mass of the caisson.
        value between 0 and 1. Required argument for RC and CC.
    rho_c : float
        density of concrete [kg/m³], required argument for RC and CC.
    rho_fill : float
        density of the fill material [kg/m³], required argument for RC
        and CC.
    Bm : float
        width of the berm [m], required argument for RC and CC.
    hb : float
        height of the foundation layer [m], required argument for RC
        and CC.
    mu : float
        friction factor between the caisson and the foundation [-],
        required argument for RC and CC.
    BermMaterial : obj
        should be a :py:class:`RockGrading` or armour unit class which
        inherits from :py:class:`ConcreteArmour`, for instance
        :py:class:`Xbloc` or :py:class:`XblocPlus`. Required argument
        for CC
    beta : float, optional, default: 0
        angle between direction of wave approach and a line normal to
        the breakwater (degrees). Optional argument for RC and CC,
        default value is 0.
    slope_toe : tuple, optional, default: (2,3)
        Slope of the armour layer (V,H). For example a slope of 2V:3H
        is defined as (2,3). Optional argument for RRM and CRM, default
        value is (2,3)
    B_toe : float, optional, default: None
        width of the top of the toe in meters. By default the width of
        toe is taken as 3 * Dn50_toe. Optional argument for RRM and
        CRM.
    vdm : {min, max, avg}, optional, default: max
        value to return in case both the deep and shallow water formula
        are valid. min for the lowest value, max for the highest value
        and avg for the average value, default is max. Optional argument
        for RRM.
    layers_rock : int, optional, default: 2
        number of layer in the armour layer for a rubble mound breakwater
        with rock as armour layer. Optional argument for RRM and RC,
        default value is 2
    layers_units : int, optional, default: 1
        number of layer in the armour layer for a rubble mound
        breakwater with concrete armour units as armour layer. Optional
        argument for RRM and CC, default value is 1
    layers_underlayer : int, optional, default: 2
        number of layers in the underlayer. Optional argument for RRM
        and CRM
    filter_rule : {'Rock', 'Xbloc', 'XblocPlus'}, optional, default: None
        filter rule to use for the substructure of the breakwater, for
        Rock, Xbloc and XblocPlus the correct filter rule is
        automatically selected. In case another type of armour layer is
        used one of these filter rules must be chosen. Optional argument
        for CRM, required if armour unit is not Xbloc or XblocPlus,
        default value is None
    pe_max : float, optional, default: 500
        maximum value of the bearing pressure at the heel of the caisson.
        Default value is set to 500 kPa, Goda (2000) advises a value
        between 400 and 500 kPa. Optional parameter for RC and CC
    SF_sliding : float, optional, default: 1.2
        safety factor against sliding. Default value according to
        Goda (2000). Optional argument for RC and CC, default value is 1.2
    SF_turning : float, optional, default: 1.2
        safety factor against sliding. Default value according to
        Goda (2000). Optional argument for RC and CC, default value is 1.2
    slope_foundation : tuple, optional, default: (2,3)
        Slope of the armour layer (V,H). For example a slope of 2V:3H
        is defined as (2,3). Optional argument for RC and CC, default
        is (2,3)
    lambda_ : list, optional, default: [1,1,1]
        modification factors of Takahasi (2002) for alternative
        monolithic breakwater. Input must be
        \\lambda_= [:math:`\\lambda_1, \\lambda_2, \\lambda_3`].
        Optional argument for RC and CC, default value is [1,1,1]
    Soil : :py:class:`Soil`, optional, default: None
        by default Soil is None, which means that the geotechnical checks
        are not performed. By specifying a Soil object, the geotechnical
        checks are automatically performed. Optional argument for all
        structures
    phi : float, optional, default: 40
        internal friction angle of rock [degrees]. Optional argument for
        RRM and CRM.

    Attributes
    ----------
    df : pd.DataFrame
        DataFrame with all generated concepts
    """
    def __init__(
            self, structure, LimitState, rho_w, slope_foreshore, Grading,
            Soil=None, safety=1, **kwargs):
        """ See help(Configurations) for more info """
        # add LimitState, rho_w and slope_foreshore to kwargs input
        # so that input can be validated
        kwargs['LimitState'] = LimitState
        kwargs['rho_w'] = rho_w
        kwargs['slope_foreshore'] = slope_foreshore
        kwargs['safety'] = safety
        kwargs['Grading'] = Grading
        kwargs['Soil'] = Soil

        # set grading as attribute for adding cost
        self._Grading = Grading

        # convert the input of structure to a list
        if isinstance(structure, list):
            # must be a list so no change
            structure = structure
        elif isinstance(structure, str):
            # convert single input to list
            structure = [structure]

        # set empty configs and vkwargs
        RM_config = {}
        RM_vkwargs = {}
        C_config = {}
        C_vkwargs = {}

        # check if the input is correct for the type of structure
        # and set default values if a default value must be set
        if 'RRM' in structure and 'CRM' in structure:
            RRM_compute = True
            CRM_compute = True
            CRMR_compute = False

            RM_vkwargs = _RM_vkwargs(type='both')
            RM_config = _process_kwargs(kwargs=kwargs, vkwargs=RM_vkwargs)

            # unpack constant arguments for both
            beta = RM_config['beta']
            phi = RM_config['phi']

            # unpack constant arguments for RRM
            N = RM_config['N']
            layers_rock = RM_config['layers_rock']
            vdm = RM_config['vdm']

            # unpack constant arguments for CRM
            ArmourUnit = RM_config['ArmourUnit']
            layers_units = RM_config['layers_units']
            filter_rule = RM_config['filter_rule']

        elif 'RRM' in structure:
            RRM_compute = True
            CRM_compute = False
            CRMR_compute = False

            RM_vkwargs = _RM_vkwargs(type='Rock')
            RM_config = _process_kwargs(kwargs=kwargs, vkwargs=RM_vkwargs)

            # unpack constant arguments for RRM
            beta = RM_config['beta']
            N = RM_config['N']
            layers_rock = RM_config['layers_rock']
            vdm = RM_config['vdm']
            phi = RM_config['phi']

        elif 'CRM' in structure:
            RRM_compute = False
            CRM_compute = True
            CRMR_compute = False

            RM_vkwargs = _RM_vkwargs(type='ArmourUnit')
            RM_config = _process_kwargs(kwargs=kwargs, vkwargs=RM_vkwargs)

            # unpack constant arguments for CRM
            beta = RM_config['beta']
            ArmourUnit = RM_config['ArmourUnit']
            layers_units = RM_config['layers_units']
            filter_rule = RM_config['filter_rule']
            phi = RM_config['phi']

        elif 'CRMR' in structure:
            RRM_compute = False
            CRM_compute = False
            CRMR_compute = True

            RM_vkwargs = _RM_vkwargs(type='ArmourUnit')
            RM_config = _process_kwargs(kwargs=kwargs, vkwargs=RM_vkwargs)

            # unpack constant arguments for CRM
            beta = RM_config['beta']
            ArmourUnit = RM_config['ArmourUnit']
            layers_units = RM_config['layers_units']
            filter_rule = RM_config['filter_rule']
            phi = RM_config['phi']

        elif 'RRM' not in structure and 'CRM' not in structure and 'CRMR' not in structure:
            RRM_compute = False
            CRM_compute = False
            CRMR_compute = False

        # input for caisson in structure
        if 'RC' in structure and 'CC' in structure:
            RC_compute = True
            CC_compute = True

            C_vkwargs = _C_vkwargs(type='both')
            C_config = _process_kwargs(kwargs=kwargs, vkwargs=C_vkwargs)

            # unpack constant arguments for both
            rho_c = C_config['rho_c']
            rho_fill = C_config['rho_fill']
            mu = C_config['mu']
            SF_sliding = C_config['SF_sliding']
            SF_turning = C_config['SF_turning']
            beta = C_config['beta']
            lambda_ = C_config['lambda_']
            filter_rule = C_config['filter_rule']
            pe_max = C_config['pe_max']

            # unpack constant arguments for RC
            layers_rock = C_config['layers_rock']

            # unpack constant arguments for CC
            layers_units = C_config['layers_units']
            BermMaterial = C_config['BermMaterial']

        elif 'RC' in structure:
            RC_compute = True
            CC_compute = False

            C_vkwargs = _C_vkwargs(type='Rock')
            C_config = _process_kwargs(kwargs=kwargs, vkwargs=C_vkwargs)

            # unpack constant arguments for both
            rho_c = C_config['rho_c']
            rho_fill = C_config['rho_fill']
            mu = C_config['mu']
            SF_sliding = C_config['SF_sliding']
            SF_turning = C_config['SF_turning']
            beta = C_config['beta']
            lambda_ = C_config['lambda_']
            filter_rule = C_config['filter_rule']
            pe_max = C_config['pe_max']

            # unpack constant arguments for RC
            layers_rock = C_config['layers_rock']

        elif 'CC' in structure:
            RC_compute = False
            CC_compute = True

            C_vkwargs = _C_vkwargs(type='ArmourUnit')
            C_config = _process_kwargs(kwargs=kwargs, vkwargs=C_vkwargs)

            # unpack constant arguments for both
            rho_c = C_config['rho_c']
            rho_fill = C_config['rho_fill']
            mu = C_config['mu']
            SF_sliding = C_config['SF_sliding']
            SF_turning = C_config['SF_turning']
            beta = C_config['beta']
            lambda_ = C_config['lambda_']
            filter_rule = C_config['filter_rule']
            pe_max = C_config['pe_max']

            # unpack constant arguments for CC
            layers_units = C_config['layers_units']
            BermMaterial = C_config['BermMaterial']

        elif 'RC' not in structure and 'CC' not in structure:
            RC_compute = False
            CC_compute = False

        # check if a computation must be made
        if (not RRM_compute and not CRM_compute and not CRMR_compute
                and not RC_compute and not CC_compute):
            # noting to compute, raise error
            raise NotSupportedError(
                (f'{structure} has not been implemented, supported structures'
                  ' are \'RRM\', \'CRM\', \'RC\' and \'CC\''))

        # initialize a df for to store the _get_concept_set
        self.df = pd.DataFrame()

        # design all concepts for a rubble mound breakwater
        # get all possible combinations of the varying arguments
        RM_varying, RM_num_combinations = self._get_combinations(
            vkwargs=RM_vkwargs, config=RM_config)

        for i in range(RM_num_combinations):
            # in first iteration check if computation is needed
            # and unpack arguments if needed
            if i == 0:
                if not RRM_compute and not CRM_compute and not CRMR_compute:
                    # break loop if computation is not needed
                    break

                if RRM_compute and CRM_compute:
                    # both thus double the number of computations
                    num = 2*RM_num_combinations
                else:
                    # only design 1 structure
                    num = RM_num_combinations

                # initialize progress bar for rubble mound computations
                RM_bar = ProgressBar(
                    number=num, task='Computing Rubble Mound')

            # set id of bw
            id = i + 1

            # get the current concept
            concept = self._get_concept_set(configs=RM_varying, id=id)

            # unpack varying arguments same for RRM and CRM
            B = concept['B']
            Dn50_core = concept['Dn50_core']
            B_toe = concept['B_toe']
            slope = concept['slope']
            slope_toe = concept['slope_toe']

            if RRM_compute:

                try:
                    with catch_warnings(record=True) as w:
                        RM_rock = RockRubbleMound(
                            slope=slope, slope_foreshore=slope_foreshore,
                            rho_w=rho_w, B=B, N=N, LimitState=LimitState,
                            Grading=Grading, Dn50_core=Dn50_core,
                            safety=safety, slope_toe=slope_toe, B_toe=B_toe,
                            layers=layers_rock, vdm=vdm, Soil=Soil, phi=phi,
                            id=id)

                except RockGradingError:
                    # if there is no rock class in the grading that
                    # can satisfy the computed Dn50 of the armour layer
                    RM_rock = None

                # save the concept to a temporary df and append to df
                temp_df = pd.DataFrame(data={'type': ['RRM'],
                                             'id': [id],
                                             'concept': [RM_rock],
                                             'B': [B],
                                             'Dn50_core': [Dn50_core],
                                             'B_toe': [B_toe],
                                             'slope': [slope],
                                             'slope_toe': [slope_toe],
                                             'warnings': [w]})
                self.df = self.df.append(temp_df, ignore_index=True, sort=True)

                RM_bar.next()

            if CRM_compute:
                try:
                    with catch_warnings(record=True) as w:
                        RM_units = ConcreteRubbleMound(
                            slope=slope, slope_foreshore=slope_foreshore, B=B,
                            rho_w=rho_w, LimitState=LimitState, safety=safety,
                            Grading=Grading, ArmourUnit=ArmourUnit, phi=phi,
                            Dn50_core=Dn50_core, slope_toe=slope_toe,
                            B_toe=B_toe, layers=layers_units, Soil=Soil, id=id,
                            filter_rule=filter_rule)

                except ArmourUnitsError:
                    # if there is no class of armour unit that
                    # can satisfy the computed Dn50 of the armour layer
                    units_bw = None

                # save the concept to a temporary df and append to df
                temp_df = pd.DataFrame(data={'type': ['CRM'],
                                             'id': [id],
                                             'concept': [RM_units],
                                             'B': [B],
                                             'Dn50_core': [Dn50_core],
                                             'B_toe': [B_toe],
                                             'slope': [slope],
                                             'slope_toe': [slope_toe],
                                             'warnings': [w]})
                self.df = self.df.append(temp_df, ignore_index=True, sort=True)

                RM_bar.next()

            if CRMR_compute:
                try:
                    with catch_warnings(record=True) as w:
                        RM_units = ConcreteRubbleMoundRevetment(
                            slope=slope, slope_foreshore=slope_foreshore, B=B,
                            rho_w=rho_w, LimitState=LimitState, safety=safety,
                            Grading=Grading, ArmourUnit=ArmourUnit, phi=phi,
                            Dn50_core=Dn50_core, slope_toe=slope_toe,
                            B_toe=B_toe, layers=layers_units, Soil=Soil, id=id,
                            filter_rule=filter_rule)

                except ArmourUnitsError:
                    # if there is no class of armour unit that
                    # can satisfy the computed Dn50 of the armour layer
                    units_bw = None

                # save the concept to a temporary df and append to df
                temp_df = pd.DataFrame(data={'type': ['CRMR'],
                                             'id': [id],
                                             'concept': [RM_units],
                                             'B': [B],
                                             'Dn50_core': [Dn50_core],
                                             'B_toe': [B_toe],
                                             'slope': [slope],
                                             'slope_toe': [slope_toe],
                                             'warnings': [w]})
                self.df = self.df.append(temp_df, ignore_index=True, sort=True)

                RM_bar.next()
            if id == RM_num_combinations:
                RM_bar.finish()

        # design all concepts for a caisson breakwater
        # get all possible combinations of the varying arguments
        C_varying, C_num_combinations = self._get_combinations(
            vkwargs=C_vkwargs, config=C_config)

        for i in range(C_num_combinations):
            # in first iteration check if computation is needed
            # and unpack arguments if needed
            if i == 0:
                if not RC_compute and not CC_compute:
                    # break loop if computation is not needed
                    break

                if RC_compute and CC_compute:
                    # both thus double the number of computations
                    num = 2*C_num_combinations
                else:
                    # only design 1 structure
                    num = C_num_combinations

                if RRM_compute or CRM_compute:
                    # set task name for nice allignment with RM
                    task_name = 'Computing Caisson     '
                else:
                    # normal task name
                    task_name = 'Computing Caisson'

                # initialize progress bar for caisson computations
                C_bar = ProgressBar(number=num, task=task_name)

            # set id of bw
            id = i + 1

            # get the current concept
            concept = self._get_concept_set(configs=C_varying, id=id)

            # unpack varying arguments
            Pc = concept['Pc']
            Bm = concept['Bm']
            hb = concept['hb']
            slope_foundation = concept['slope_foundation']

            if RC_compute:
                # design with rock as armour layer for the foundation
                try:
                    with catch_warnings(record=True) as w:
                        C_rock = Caisson(
                            Pc=Pc, rho_c=rho_c, rho_fill=rho_fill, beta=beta,
                            rho_w=rho_w, Bm=Bm, hb=hb, layers=layers_rock,
                            BermMaterial=Grading, mu=mu, LimitState=LimitState,
                            safety=safety, slope_foreshore=slope_foreshore,
                            SF_sliding=SF_sliding, SF_turning=SF_turning,
                            slope_foundation=slope_foundation, lambda_=lambda_,
                            filter_rule=filter_rule, Grading=Grading, id=id,
                            pe_max=pe_max, Soil=Soil)

                except RockGradingError:
                    # if there is no rock class in the grading that
                    # can satisfy the computed Dn50 of the armour layer
                    C_rock = None

                # save the concept to a temporary df and append to df
                temp_df = pd.DataFrame(data={'type': ['RC'],
                                             'id': [id],
                                             'concept': [C_rock],
                                             'Pc': [Pc],
                                             'Bm': [Bm],
                                             'hb': [hb],
                                             'slope_foundation': [slope_foundation],
                                             'warnings': [w]})
                self.df = self.df.append(temp_df, ignore_index=True, sort=True)

                C_bar.next()

            if CC_compute:

                try:
                    with catch_warnings(record=True) as w:
                        C_units = Caisson(
                            Pc=Pc, rho_c=rho_c, rho_fill=rho_fill, rho_w=rho_w,
                            Bm=Bm, hb=hb, layers=layers_units, mu=mu, beta=beta,
                            BermMaterial=BermMaterial, LimitState=LimitState,
                            slope_foreshore=slope_foreshore, safety=safety,
                            SF_sliding=SF_sliding, SF_turning=SF_turning,
                            slope_foundation=slope_foundation, lambda_=lambda_,
                            filter_rule=filter_rule, Grading=Grading, id=id,
                            pe_max=pe_max, Soil=Soil)

                except ArmourUnitsError:
                    # if there is no class of armour unit that
                    # can satisfy the computed Dn50 of the armour layer
                    C_units = None

                # save the concept to a temporary df and append to df
                temp_df = pd.DataFrame(data={'type': ['CC'],
                                             'id': [id],
                                             'concept': [C_units],
                                             'Pc': [Pc],
                                             'Bm': [Bm],
                                             'hb': [hb],
                                             'warnings': [w]})
                self.df = self.df.append(temp_df, ignore_index=True, sort=True)

                C_bar.next()

            if id == C_num_combinations:
                C_bar.finish()

    @staticmethod
    def _get_concept_set(configs, id):
        """ get unique set of parameters and values """
        index = id - 1

        unique_set = {}

        for param, val in configs.items():
            unique_set[param] = val[index]

        return unique_set

    @staticmethod
    def _get_combinations(vkwargs, config):
        # get the varying and constant parameters
        varying = {}
        fixed = {}

        for param, val in config.items():
            if vkwargs[param]['Constant']:
                # constant parameter, already in RM_config
                continue
            else:
                # varying parameter
                if isinstance(val, np.ndarray):
                    # also given as varying parameter
                    varying[param] = val
                else:
                    # not set as varying thus constant
                    fixed[param] = val

        # create all possible combinations of the varying parameters
        # get the parameters and number of parameters
        parameters = list(varying.keys())
        num_parameters = len(parameters)

        # empty dict to store the parameter set of each concept
        configs = {}

        # create value_index with the index of the value in the list of the dict
        # value_index with 0 is the first concept
        value_index = [0] * num_parameters
        first_concept = True
        current_combination = 0

        while True:
            # check for exit condition
            # first compute maximum number of combinations
            if first_concept:
                max_combinations = 1
                for param, values in varying.items():
                    combinations = len(values)
                    max_combinations = max_combinations * combinations

            if current_combination == max_combinations:
                break

            # saving the concept in configs
            for i in range(num_parameters):
                current_parameter = parameters[i]
                if first_concept:
                    # first time a list for each key must be made
                    configs[current_parameter] = []
                # save the value of the current key
                value = varying[current_parameter][value_index[i]]
                configs[current_parameter].append(value)

            # first concept has been generated so set to False
            first_concept = False


            # set loop variables, current_column_index starts from the right
            change = True
            current_column_index = num_parameters - 1

            while change and current_column_index >= 0:
                # get the length of the current parameter (how many values)
                current_parameter = parameters[current_column_index]
                max_values = len(varying[current_parameter])

                # check if all variants for this level have been made
                if (value_index[current_column_index] + 1) > max_values-1:
                    # all variants for the level have been added
                    # so set the last index to zero
                    value_index[current_column_index] = 0

                    # Change the upper variable by one
                    # We need to increment the immediate upper level loop by one
                    change = True
                else:
                    # add one to the index for the next combination
                    value_index[current_column_index] += 1

                    # set the change to False so that the loop stops for the
                    # current level
                    change = False

                # move one column to the left
                current_column_index -= 1

            current_combination += 1

        # check if there are fixed values
        if fixed:
            # add the fixed parameters to the configs
            for param, val in fixed.items():
                values = [val for i in range(max_combinations)]
                configs[param] = values

        # return the combinations and number of combinations
        return configs, max_combinations

    def add_cost(
            self, type = 'Material', equipment= None, core_price=None, unit_price=None, concrete_price=None,
            fill_price=None, transport_cost=None, investment=None,
            length=None):
        """ Compute the cost of each concept either CO2 or material cost

        Compute the cost of each concept and add the cost to the
        :py:attr:`df`. The cost of the rocks must be specified in the
        RockGrading. If transport cost are not included in the price of
        rocks or core_price it can be given with the argument
        transport_cost. For a Caisson breakwater it is possible to
        specify the investment for renting a dry dock, the investment
        is divided through the length of the breakwater to get the
        investment cost per meter.

        .. note::
           The transport_cost are not added to the price of the armour
           layer. The assumption has been made that the cost of
           producing and transporting the armour units is included in
           the unit_price.

        Parameters
        ----------
        type: {'Material', 'CO2'}
            which cost type is added
        Equipment: lst
            With equipment from Equipment Class
        core_price : float, optional, default: None
            cost of the core material per m³, required for RRM and CRM
        unit_price : float, optional, default: None
            the cost of an armour unit per m³, required for CRM and CC
        concrete_price : float, optional, default: None
            price of concrete per m³, required for RC and CC
        fill_price : float, optional, default: None
            price of the fill material per m³, required for RC and CC
        transport_cost : float, optional, default: None
            the cost to transport a m³ of rock from the quarry to the
            project location
        investment : float
            the investment required to rent a dry dock
        length : float
            length of the breakwater [m]
        """
        # make dict of the cost for validation
        cost = {
            'core_price': core_price,
            'unit_price': unit_price,
            'concrete_price': concrete_price,
            'fill_price': fill_price}

        # check if all required cost have been given
        for structure in self.df.type.unique():
            # validate cost
            _process_cost(structure, type, cost, self._Grading)

        # set list to store cost in
        computed_cost = []
        # iterate over the generated concepts
        for i, row in self.df.iterrows():
            # check if concept is not None
            if row.concept is None:
                # not a valid concept
                computed_cost.append(None) #Wat gebeurt hier. Dit zorgt voor een error bij cost_influence
            else:
                # valid concept
                # check types and compute price
                if row.type == 'RRM':
                    price = row.concept.cost(
                        *row.concept.variantIDs, type = type, equipment = equipment, core_price= core_price,
                        transport_cost=transport_cost, output='variant')

                elif row.type == 'CRM':
                    price = row.concept.cost(
                        *row.concept.variantIDs, type = type, equipment = equipment, core_price= core_price,
                        unit_price=unit_price, transport_cost=transport_cost,
                        output='variant')

                elif row.type == 'CRMR':
                    price = row.concept.cost(
                        *row.concept.variantIDs, type = type, equipment = equipment, core_price= core_price,
                        unit_price=unit_price, transport_cost=transport_cost,
                        output='variant')

                elif row.type == 'RC' or row.type == 'CC':
                    # check if investment cost must be added
                    if investment is not None and length is not None:
                        # add investment cost
                        row.concept.dry_dock(investment, length)

                    price = row.concept.cost(
                        *row.concept.variantIDs, type = type, equipment = equipment,core_price= core_price, concrete_price=concrete_price,
                        fill_price=fill_price, unit_price=unit_price)

                else:
                    raise NotSupportedError(f'{row.type} is not supported')

                # add cost to list
                computed_cost.append(price)

        # add column to the df for either material or CO2
        if type == 'Material':
            self.df['material_cost'] = computed_cost
        if type == 'CO2':
            self.df['CO2_cost'] = computed_cost

    def to_design_explorer(
            self, params, mkdir='DesignExplorer', slopes='angles',
            merge_Bm=True, merge_slope_toe=True):
        """ Export concepts to Design Explorer 2

        Creates a folder that can be used in Design Explorer 2, the
        folder consists of the cross sections of all concepts and a csv
        file with the data of the concepts. Parameters supported for
        export can be seen in table 3. To use the folder in Design
        Explorer 2, follow these steps:

        - Upload the folder to your Google Drive
        - Share the folder and get the shareable link
        - Go to http://tt-acm.github.io/DesignExplorer/
        - Click on *Get Data* and paste the shareable link below: *From
          the cloud*
        - Click on *Load Data*
        - Enjoy exploring!

        Table 3: possible parameter to export

        +--------------------------+------------+------------+
        | Parameter                | RRM + CRM  |  RC + CC   |
        +==========================+============+============+
        | material_cost            |     o      |     o      |
        +----------------------------------------------------+
        | CO2_cost                 |     o      |     o      |
        +--------------------------+------------+------------+
        | B                        |     o      |     o      |
        +--------------------------+------------+------------+
        | Rc                       |     o      |     o      |
        +--------------------------+------------+------------+
        | computed Dn50 armour     |     o      |     o      |
        +--------------------------+------------+------------+
        | class armour             |     o      |     o      |
        +--------------------------+------------+------------+
        | class Dn50 armour        |     o      |     o      |
        +--------------------------+------------+------------+
        | computed Dn50 underlayer |     o      |     o      |
        +--------------------------+------------+------------+
        | class underlayer         |     o      |     o      |
        +--------------------------+------------+------------+
        | class Dn50 underlayer    |     o      |     o      |
        +--------------------------+------------+------------+
        | slope                    |     o      |            |
        +--------------------------+------------+------------+
        | slope_toe                |     o      | o :sup:`1` |
        +--------------------------+------------+------------+
        | Dn50_core                |     o      |            |
        +--------------------------+------------+------------+
        | B_toe                    |     o      |            |
        +--------------------------+------------+------------+
        | computed Dn50 filter     |     o      |            |
        +--------------------------+------------+------------+
        | class filter             |     o      |            |
        +--------------------------+------------+------------+
        | class Dn50 filter        |     o      |            |
        +--------------------------+------------+------------+
        | Pc                       |            |     o      |
        +--------------------------+------------+------------+
        | hb                       |            |     o      |
        +--------------------------+------------+------------+
        | h_acc                    |            |     o      |
        +--------------------------+------------+------------+
        | Bm                       | o :sup:`2` |     o      |
        +--------------------------+------------+------------+
        | slope_foundation         |            |     o      |
        +--------------------------+------------+------------+
        | UC :sup:`3`              |            |     o      |
        +--------------------------+------------+------------+

        | :sup:`1` slope_foundation is interpreted as slope_toe if
          merge_slope_toe is set to True
        | :sup:`2` B_toe is interpreted as Bm if merge_Bm is set to True
        | :sup:`3` UC is the unity check for the bearing capacity of the
          subsoil

        Parameters
        ----------
        params : list
            list of parameters for the design explorer, parameter must
            be str. See table 3 for the parameters that can be exported/
        mkdir : str, optional, default: concepts
            creates a folder in which all cross sections are saved,
            upload the created folder to your Google Drive to use
            Design Explorer 2
        slopes : {angles, tuples}, optional, default: angles
            how the slopes must be exported. tuples will export the
            slope as (V,H) and angles as an angle in degrees
        merge_Bm : bool, optional, default: True
            True if Bm must be merged with B_toe when a Rubble Mound and
            Vertical breakwater have been designed, False if you do not
            want to merge the columns.
        merge_slope_toe : bool, optional, default: True
            True if slope_foundation must be merged with slope_toe when
            a Rubble Mound and Vertical breakwater have been designed,
            False if you do not want to merge the columns.

        Raises
        ------
        KeyError
            if the cost are asked to export but not yet specified
        """
        # set empty df for export
        to_export = pd.DataFrame()

        # check if RRM and CRM are in structure
        designed_structures = self.df.type.unique()
        if 'CRM' and 'RRM' and 'CRMR' in designed_structures:
            format_class_as_string = True
        else:
            format_class_as_string = False

        # set dir int variable for generating a new dir
        dir_int = 0

        # while loop to create dir with unique name
        while True:
            # set name of the new_dir
            if dir_int == 0:
                # first iteration, so use specified name
                new_dir = mkdir
            else:
                # specified name already exists so add int to name of dir
                new_dir = f'{mkdir} {dir_int}'

            # check if dir already exists
            if os.path.exists(new_dir):
                # dir already exists, check if it has files
                if len(os.listdir(new_dir)) == 0:
                    # no files thus data can be added, break loop
                    break
                else:
                    # increase dir_int with 1 to make new dir
                    dir_int += 1
            else:
                # dir does not exists, so make one and break
                os.mkdir(new_dir)
                break

        # check if the name of dir is different from the specified
        if mkdir is not new_dir:
            # name different, print user warning to notify user
            user_warning(
                (f'directory {mkdir} already exists, therefore a new directory'
                 f' has been made. The files can be found in {new_dir}'))

        # set up progress bar, and CaseNo
        num_concepts = self.df.shape[0]
        bar = ProgressBar(number=num_concepts, task='Saving')
        CaseNo = 1

        # add list to store all CaseNo in
        all_CaseNo = []

        # start the export
        for index, row in self.df.iterrows():
            # check if concept exists
            if row.concept is None: #Als het concept None is dan wordt het eruit gehaald
                # go to next concept
                bar.next()
                all_CaseNo.append([CaseNo])
                continue

            # compute number of variants for this concept
            num_concepts = len(row.concept.variantIDs)

            # store CaseNo of current concept in a list
            temp_CaseNo = []

            for id in row.concept.variantIDs:
                # save name of the cross section
                file_name = f'{row.type}.{row.id}'

                if num_concepts > 1:
                    # add id if more than 1 concept
                    file_name = f'{file_name}{id}'

                save_name = f'{new_dir}/{file_name}'

                # save cross section of the current concept
                row.concept.plot(id, save_name=save_name)

                data = {'CaseNo': [CaseNo], 'type': [row.type]}

                # get the values from the variant and add to data
                variant = row.concept.get_variant(variantID=id)
                to_explorer = _DE_params(
                    args=params, variant=variant, row=row, concept=row.concept,
                    structure=row.type, slopes=slopes,
                    change_CRM_class=format_class_as_string)
                data.update(to_explorer)

                # check if material cost must be included
                if 'material_cost' in params:
                    # check if cost have been added
                    if 'material_cost' in self.df.columns:
                        # add cost to data
                        data['material_cost'] = row.material_cost[id]
                    else:
                        raise KeyError(
                            'Material cost have not been added, use add_cost to add '
                            'cost to the df')

                # check if CO2 cost must be included
                if 'CO2_cost' in params:
                    # check if cost have been added
                    if 'CO2_cost' in self.df.columns:
                        # add cost to data
                        data['CO2_cost'] = row.CO2_cost[id]
                    else:
                        raise KeyError(
                            'CO2 cost have not been added, use add_cost to add '
                            'cost to the df')

                # add image to data
                data['img'] = [f'{file_name}.png']

                # create a df of data and add to the df to_export
                export_row = pd.DataFrame(data=data)
                to_export = to_export.append(
                    export_row, ignore_index=True, sort=True)

                # add CaseNo to temp CaseNo and increase with 1
                temp_CaseNo.append(CaseNo)
                CaseNo += 1

            # add all stored CaseNo of the current concept in a list
            all_CaseNo.append(temp_CaseNo)

            bar.next()

        bar.finish()

        # add all CaseNo to the df so that concept can be selected by CaseNo
        self.df['CaseNo'] = all_CaseNo

        # post process for to_export
        # if CRM and RRM for class armour values must be sorted
        if format_class_as_string:
            to_export.sort_values(by=['type'], inplace=True)

        # check if all given params are in the columns
        columns = list(to_export)
        not_used_params = np.setdiff1d(
            params, columns, assume_unique=True).tolist()
        if any(not_used_params):
            # remove parameters not used from params
            column_names = [x for x in params if x not in not_used_params]

            # raise UserWarning that not all given params have been exported
            skipped = ', '.join(not_used_params)
            user_warning(
                f'Not all params have been exported, skipped: {skipped}')
        else:
            # all given params have been used
            column_names = params

        # add CaseNo and type as the left most column_names
        # and make img the right most column_names
        column_names[0:0] = ['CaseNo', 'type']
        column_names.insert(len(column_names), 'img')

        # restructure df with order of column_names
        to_export = to_export.loc[:, column_names]

        # merge B_toe with Bm if set to True
        if 'B_toe' in column_names and 'Bm' in column_names and merge_Bm:
            to_export.Bm.fillna(to_export.B_toe, inplace=True)
            del to_export['B_toe']

        # merge slope_foundation with slope_toe if set to True
        if ('slope_toe' in column_names and 'slope_foundation' in column_names
                and merge_slope_toe):
            to_export.slope_toe.fillna(
                to_export.slope_foundation, inplace=True)
            del to_export['slope_foundation']

        # save to_export df to a excel file
        excel_save_name = f'{new_dir}/data.csv'
        to_export.to_csv(excel_save_name, index=False)

        print(f'folder {new_dir} is ready for Design Explorer 2')

    def get_concept(self, id=None, CaseNo=None):
        """ Get the specified concept

        Parameters
        ----------
        id : str, optional, default: None
            id of the concept, the id of a concept for a rubble mound
            breakwater out of rock is for instance: RRM.1
        CaseNo : int, optional, default: None
            CaseNo if the concept, only available if an export for the
            Design Explorer has been made with
            :py:meth:`to_design_explorer`

        Returns
        -------
        concept : obj
            a breakwater concept, for instance :py:obj:`RockRubbleMound`

        Raises
        ------
        InputError
            if a concept is not selected with CaseNo or id, or if there
            is no concept with the specified CaseNo or id
        KeyError
            if a concept can't be selected by CaseNo because the method
            :py:meth:`to_design_explorer` has not yet been used, and the
            CaseNo have therefore not yet been added to :py:attr:`df`
        """
        if id and CaseNo is not None:
            # cannot select concept by both id and CaseNo
            raise InputError(
                'Concept must be selected with id or CaseNo, not both')

        elif id is not None:
            # select concept by id
            type = id.split('.')[0]
            id = id.split('.')[1]

            row_ids = self.df[self.df['id'] == int(id)]
            row = row_ids[row_ids['type'] == type]

            msg = f'id {type}.{id}'

        elif CaseNo is not None:
            # check if CaseNo is a column in df,
            # is only added if the method to_design_explorer is used
            if 'CaseNo' not in self.df.columns:
                raise KeyError(
                    ('CaseNo will only be added to the df if an export for '
                     'the design explorer has been made, run '
                     'to_design_explorer or select by id'))

            # get the row of the df where the given CaseNo is in column CaseNo
            row = self.df[self.df.apply(
                lambda x: CaseNo in x['CaseNo'], axis=1)]

            msg = f'CaseNo {CaseNo}'

        else:
            # CaseNo and id is not given thus raise error
            raise InputError(
                ('No concept could be selected as neither an id nor CaseNo '
                 ' has been given'))

        # check if row has values
        if row.empty:
            # df is empty, thus CaseNo was invalid
            raise InputError(f'There is no concept with {msg}')
        else:
            # return the concept
            return row['concept'].values[0]

    def to_breakwaters(self, save_name):
        """ Save the df with designs in a pickle object

        Method will save the df with all breakwater concepts in a pickle
        object with extension .breakwaters.

        Parameters
        ----------
        save_name : str
            name of the file
        """
        # save to pickle
        file_name = f'{save_name}.breakwaters'

        file = open(file_name, 'wb')
        pickle.dump(self, file)

        print(f'saved to {file_name}')

    def show_warnings(self):
        """ Print all warnings encountered during the design

        Method prints a table containing the unique warning messages,
        with the time the warning was encountered during the design.
        """
        all_warnings = self.df.warnings

        unique_warnings = {}

        # iterate over de column in the df
        for warnings in all_warnings:
            # iterate over the list of warnings from each concept
            for warning in warnings:
                msg = str(warning.message)
                # check if warning is already encountered
                if msg in unique_warnings:
                    # not unique, thus add one to the count
                    unique_warnings[msg] += 1
                else:
                    # warning is unique, thus add to unique-wa
                    unique_warnings[msg] = 1

        # print the table
        table = [[key, count] for key, count in unique_warnings.items()]
        print(tabulate(table, headers=['Warning', 'Count'], tablefmt='github'))

    def cost_influence(self, type_):
        """ Plot the influence of the varying parameters on the cost

        Method to see the influence of the varying parameters on the
        cost of a concept. If more than 1 parameter has been specified
        as a varying parameter the x-axis is normalised (from 0 to 1)
        so that all parameters can be plotted on the same axis.

        type_: {'Material, 'CO2'}
            which cost influence is analysed

        .. note::
           When more than one varying parameter is specified several
           lines must be plotted. To show the influence of one varying
           parameter on the cost the other varying parameters are set
           to their lower bound value.
        """
        # headers to exclude
        exclude = ['concept', 'id', 'type', 'warnings']

        cost_var = None

        if type_ == 'Material':
            cost_var = 'material_cost'
        if type_ == 'CO2':
            cost_var = 'CO2_cost'

        # make dict to the line data for plotting
        lines = {}

        for structure in self.df.type.unique():
            # make df where the type equals the current structure
            df = self.df[self.df.type == structure]
            df = df.drop(columns=exclude)
            # iterate over the columns of the df
            for parameter in df.columns:
                # check if parameter is slope
                if 'slope' in parameter:
                    # check if not all nan
                    if not df[parameter].isnull().all():
                        # convert tuples to floats
                        df[parameter] = df[parameter].apply(
                            lambda x: np.round(np.arctan(x[0]/x[1])*180/np.pi, 2))

                # check if not a cost parameter or in exclude
                if parameter not in ['material_cost', 'CO2_cost'] and parameter not in exclude:
                    # get the unique values
                    unique = df[parameter].unique()
                    # check if only 1 unique value in the column
                    if len(unique) > 1:
                        # varying parameter, add to lines dict
                        lines[parameter] = {'values': [], cost_var: []}

                        # get list of the values
                        values = list(df[parameter].values)

                        # iterate over the unique values
                        for unique_val in unique:
                            # get index of value in the df, first occurance
                            index = values.index(unique_val)

                            # get the row in the df with the data
                            row = df.iloc[index]

                            # add info to dict
                            # note that the unique value is retrieved from
                            # the df and not from the list with unique
                            # values, this is done so the result can be
                            # checked, i.e. if all values are indeed unique
                            lines[parameter]['values'].append(
                                row[parameter])

                            # compute average cost and add to dict
                            #concept excluded so cost row.concept != None replaced by row.material_cost
                            if type_ == 'Material' and row.material_cost != None:
                                lines[parameter][cost_var].append(
                                    np.mean(list(row.material_cost.values())))
                            if type_ == 'CO2' and row.CO2_cost != None:
                                lines[parameter][cost_var].append(
                                    np.mean(list(row.CO2_cost.values())))

        # generate the plot

        plot = cost_influence(type= type_, lines= lines)

