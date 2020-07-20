import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

from .utils.exceptions import InputError, user_warning, NotSupportedError, RockGradingError
from .core.goda import Goda
from .core.overtopping import vertical
from .core.toe import toe_berm_stability
from .core.substructure import underlayer, _supported_armour_layers, layer_coefficient
from .core.scour import scour_protection


class Caisson:
    """ Design a (composite) vertical breakwater

    Makes a conceptual design of a vertical or composite vertical
    breakwater, with a caisson on a rubble mound foundation. The
    following computations are performed:

    - The necessary size of the armour layer of the foundation is
      designed with the modified Tanimoto formula (Takahashi, 2002).
    - The required stone size for the core of the foundation
    - The water depth in front of the caisson is computed based on the
      dimensions of the foundation and water depth
    - The crest freeboard is computed with the formulae from EurOtop
      (2018), :py:func:`vertical` is used, which automatically
      classifies the breakwater so that the correct formula is used.
    - The required width of the caisson is computed with the extended
      Goda formula (Takahasi, 2002).
    - The required width of the scour protection with Sumer and Fredsoe
      (2000). Note that a scour protection is only added if the width
      of the foundation is not sufficient.
    - If a Soil is specified the bearing capacity of the soil will also
      be checked with Brinch Hansen (1970).

    Parameters
    ----------
    Pc : float
        contribution of concrete to the total mass of the caisson.
        value between 0 and 1
    rho_c : float
        density of concrete [kg/m³]
    rho_fill : float
        density of the fill material [kg/m³]
    rho_w : float
        density of water [kg/m³]
    Bm : float
        width of the berm [m]
    hb : float
        height of the foundation layer [m]
    layers : int
        number of layers in the armour layer of the toe berm
    BermMaterial : obj
        should be a :py:class:`RockGrading` or armour unit class which
        inherits from :py:class:`ConcreteArmour`, for instance
        :py:class:`Xbloc` or :py:class:`XblocPlus`
    LimitState : :py:class:`LimitState` or list of :py:class:`LimitState`
        ULS, SLS or another limit state defined with
        :py:class:`LimitState`
    slope_foreshore : tuple
        slope of the foreshore (V, H). For example a slope of 1:100 is
        defined as (1,100)
    mu : float
        friction factor between the caisson and the foundation [-]
    safety : float, optional, default: 1
        safety factor of design (number of standard deviations from the
        mean)
    SF_sliding : float, optional, default: 1.2
        safety factor against sliding. Default value according to
        Goda (2000)
    SF_turning : float, optional, default: 1.2
        safety factor against sliding. Default value according to
        Goda (2000)
    beta : float, optional, default: 0
        angle between direction of wave approach and a line normal to
        the breakwater [degrees]
    slope_foundation : tuple, optional, default: (2,3)
        Slope of the armour layer (V, H). For example a slope of 2V:3H
        is defined as (2,3)
    lambda_ : list, optional, default: [1, 1, 1]
        modification factors of Takahasi (2002) for alternative
        monolithic breakwater. Input must be
        \\lambda_= [:math:`\\lambda_1, \\lambda_2, \\lambda_3`].
    pe_max : float, optional, default: 500
        maximum value of the bearing pressure at the heel of the caisson.
        Default value is set to 500 kPa, Goda (2000) advises a value
        between 400 and 500 kPa. [kPa]
    filter_rule : {'Rock', 'Xbloc', 'XblocPlus'}, optional, default: None
        filter rule to use for the substructure of the breakwater, for
        rock, Xbloc and XblocPlus the correct filter rule is
        automatically selected. In case another type of armour layer is
        used one of these filter rules must be chosen.
    Grading : :py:class:`RockGrading`, optional, default: None
        standard rock grading defined in the NEN-EN 13383-1 or a user
        defined rock grading. Required if the BermMaterial is not a
        :py:class:`RockGrading`.
    Soil : :py:class:`Soil`, optional, default: None
        by default Soil is None, which means that the geotechnical checks
        are not performed. By specifying a Soil object, the geotechnical
        checks are automatically performed.
    id : int, optional, default: None
        unique id of the breakwater

    Attributes
    ----------
    logger : dict
        dict of warnings and messages
    structure : dict
        dictionary with the computed Dn50, rock class, and average Dn50
        of the rock class for each layer and the toe. This dictionary
        includes all variants, use :py:meth:`get_variant` to get the
        parameters of one specific variant. Alternatively,
        :py:meth:`print_variant` can be used to print the details of
        one, multiple or all variants.
    id : int
        unique id of the breakwater
    variantIDs : list
        list with the IDs of the variants generated for this rubble
        mound breakwater.
    price : dict
        cost of each variant, generated when :py:meth:`cost` or
        :py:meth:`dry_dock` is used.
    width_scour : float
        the required length of the scour protection [m]. The distance
        is measured from the front of the caisson
    geotechnical : dict
        dictionary with the computed values from the Brinch Hansen
        equation, has keys: Fr with the resistance, N with the net
        downward force and UC, which is the unity check (N/Fr). This
        dictionary is generated when the :py:obj:`Soil` argument is
        specified.
    """

    def __init__(
            self, Pc, rho_c, rho_fill, rho_w, Bm, hb, layers, BermMaterial,
            LimitState, slope_foreshore, mu, safety=1, SF_sliding=1.2,
            SF_turning=1.2, beta=0, slope_foundation=(2,3), lambda_=[1,1,1],
            pe_max=500, filter_rule=None, Grading=None, Soil=None, id=None,
            **kwargs):
        """ See help(Caisson) for more info """
        # set logger and structure
        self.logger = {'INFO': [], 'WARNING': []}
        self.structure = {}

        # set id and variantIDs
        self.id = id
        self.variantIDs = ['a']

        # set cost Attribute
        self.price = None

        # compute relative buoyant density
        delta = (BermMaterial.rho - rho_w)/rho_w

        # convert beta from degrees to rad
        beta = beta*np.pi/180

        # compute the angle of the foreshore
        slope_foreshore = np.arctan(slope_foreshore[0]/slope_foreshore[1])

        # set top layer and filter layer
        top_layer = BermMaterial.name
        supported = _supported_armour_layers()

        if top_layer in supported:
            # supported armour layer
            filter_rule = top_layer
            if Grading is None and top_layer != 'Rock':
                # no grading for armour unit, thus raise error
                supported_armour_units = [
                    layer for layer in supported if layer is not 'Rock']
                support_out = ', '.join(supported_armour_units)
                raise InputError(
                    'Missing argument: Grading. The top layer is made out of '
                    f'{support_out}, therefore a RockGrading must be specified')
            elif top_layer == 'Rock':
                Grading = BermMaterial

        if filter_rule is None:
            supported_rules = ', '.join(supported)
            raise NotSupportedError(
                (f'Filter rule for {top_layer} is not implemented, set '
                 f'filter rule to use with filter_rule to {supported_rules}'))

        # restructure the input LimitState(s)
        self._LimitStates = []
        if isinstance(LimitState, list):
            self._LimitStates.extend(LimitState)
        else:
            self._LimitStates.append(LimitState)

        # set private attribute for plotting cross section
        self._input_arguments = {
            'Grading': Grading,
            'slope_foundation': slope_foundation,
            'armour': top_layer}

        # design the armour layer of the foundation
        # set temporary values to check changes
        Dn50_berm, state = 0, 0

        # iterate over the LimitStates
        for i, LimitState in enumerate(self._LimitStates):
            # get the wave height and submerged depth of caisson
            H13 = LimitState.get_Hs(definition='H13')
            T13 = LimitState['T13']
            h = LimitState.h
            d = h - hb

            # set values for the while loop
            Dn50 = 0
            compute_berm = True

            # compute Dn50 of the berm in a while loop,
            # because d = h - hb - layers*Dn50
            while compute_berm:
                Dn50_temp = toe_berm_stability(
                    Hs=H13, T=T13, d=d, Bm=Bm, Delta=delta, beta=beta)
                if (Dn50_temp - Dn50) > -0.05 and (Dn50_temp - Dn50) < 0.05:
                    compute_berm = False
                else:
                    Dn50 = Dn50_temp
                d = h - hb - layers*Dn50

            # check if computed Dn50 of current LimitState is larger
            # than current normative Dn50
            if Dn50 > Dn50_berm:
                # value is larger than current normative
                Dn50_berm = Dn50
                state = i

        if d <= 0:
            raise InputError(
                ('Encountered negative value for d. The chosen foundation '
                 'height (hb) is probably very close to the water depth (h)'))

        class_berm = BermMaterial.get_class(Dn50_berm)

        if top_layer == 'Rock':
            # berm material is made out of rock
            class_dn = BermMaterial.get_class_dn50(class_berm)
        else:
            # berm material is made out of concrete armour units
            class_dn = class_berm**(1/3)

        # compute d, the depth above the foundation
        d = self._LimitStates[state].h - hb - layers*Dn50_berm

        # add dimensions of armour layer to the structure
        self.structure['armour'] = {
            'computed Dn50': Dn50_berm,
            'class': class_berm,
            'class Dn50': class_dn,
            'state': state,
            'layers': layers}

        # check if additional layer is needed
        range_Dn50_u = underlayer(
            Dn_armour=class_dn, armour_layer=filter_rule,
            rho=BermMaterial.rho, rho_rock=Grading.rho)

        class_u, class_dn_u = [], []

        for dn_u in range_Dn50_u:
            rock_class = Grading.get_class(dn_u)

            if rock_class in class_u:
                # rock class is already in the underlayer so no need to
                # generate a new variant because design will be the same
                continue
            else:
                # first loop or new variant
                class_u.append(rock_class)
                class_dn_u.append(Grading.get_class_dn50(rock_class))

        self.structure['foundation'] = {'computed Dn50': range_Dn50_u,
                                        'class': class_u,
                                        'class Dn50': class_dn_u,
                                        'state': 'see armour'}

        # add message to the log if new variant was generated
        if len(class_u) == 2:
            self.logger['INFO'].append(
                'two rock classes possible for the underlayer, generated new '
                'variant b')
            self.variantIDs.append('b')

        # compute the crest height
        # set temporary values
        Rc, state_overtopping, B, m = 0, 0, 0, 0

        for i, LimitState in enumerate(self._LimitStates):
            self.logger['INFO'].append(f'computing with {LimitState.label}:')
            # get hydraulic parameters
            Hm0 = LimitState.get_Hs(definition='Hm0')
            L = LimitState.L(period='T_m_min_1', deep_water=True)
            s = LimitState.s(number='spectral')
            h = LimitState.h

            Rc_temp = vertical(
                Hm0=Hm0, q=LimitState['q'], h=h, d=d, L_m_min_1=L,
                s_m_min_1=s, safety=safety, logger=self.logger)

            # check if computed Rc of current LimitState is larger
            # than current normative Rc
            if Rc_temp > Rc:
                # if larger the normative Rc must be changed
                Rc = Rc_temp
                state_overtopping = i

            # compute/get wave heights for Goda
            H13 = LimitState.get_Hs(definition='H13')
            Hmax = LimitState['Hmax']

            # compute the height of the caisson
            h_acc = h - hb

            goda = Goda(
                Hs=H13, Hmax=Hmax, h=h, d=d, h_acc=h_acc, hc=Rc_temp, Bm=Bm,
                T=LimitState['T13'], beta=beta, rho=rho_w, logger=self.logger,
                slope_foreshore=slope_foreshore, lambda_=lambda_)

            B_temp = goda.required_width(
                Pc=Pc, rho_c=rho_c, rho_f=rho_fill, rho_w=rho_w, mu=mu,
                t=0.5, SF_sliding=SF_sliding, SF_turning=SF_turning,
                logger=self.logger)

            # check bearing capacity
            pe = goda.bearing_pressure(Pc=Pc, rho_c=rho_c, rho_fill=rho_fill)

            if pe/1000 > pe_max:
                # bearing capacity is not large enough
                B_temp = goda.bearing_pressure_width(
                    B1=B_temp, Pc=Pc, rho_c=rho_c, rho_fill=rho_fill,
                    pe_max=pe_max)

                # replace normative msg from required_width in logger
                self.logger['INFO'][-1] = ('The bearing pressure is normative'
                    ' for the computation of the width')

            m_temp = goda.mass(Pc=Pc, rho_c=rho_c, rho_fill=rho_fill)

            # check if values have changed from the current normative
            if m_temp > m:
                # set new normative values
                B = B_temp
                m = m_temp
                self.goda = goda
                state_goda = i

        # compute submerged depth of the caisson
        h_acc = self._LimitStates[state_goda].h - hb

        self.structure['caisson'] = {
            'hb': hb, 'h_acc': h_acc, 'Pc': Pc, 'd': d,
            'Rc': Rc, 'state_overtop': state_overtopping,
            'B': B, 'state_goda': state_goda, 'Bm': Bm}

        # Compute required scour protection
        self.width_scour = 0
        for LimitState in self._LimitStates:
            w = scour_protection(L=LimitState.L(period='Tm'))

            # check if larger than previous value
            if w >= self.width_scour:
                # set w as new width scour
                self.width_scour = w

        # make geotechnical design
        if Soil is not None:
            # compute the effective width of the foundation
            B_eff_sub = self._effective_width_foundation(rho_w, hb, B, m)

            # compute the forces of the foundation on the subsoil
            # compute horizontal stress per meter
            t = (self.goda.P()/(h_acc + Rc))/B_eff_sub/1000

            # compute vertical stress per meter
            p = (self.goda._dFv(m) + self._Fsill(rho_w, hb, B))/B_eff_sub/1000

            # compute the bearing capacity of the subsoil
            p_cap = Soil.brinch_hansen(
                p=p, t=t, B=B_eff_sub, L=None, q=0, rho_w=rho_w)

            # check if capacity is enough
            if p >= p_cap:
                # force is larger
                user_warning(
                    ('Bearing capacity of the soil is smaller than the '
                     'exerted stress'))

            # make geotechnical attribute and add values
            self.geotechnical = {'p_cap': p_cap, 'p': p, 'UC': p/p_cap}

    def _Fsill(self, rho_w, hb, B):
        """ Compute the downward force of the foundation """
        # slope
        V, H = 1, 1

        # get the grading
        Grading = self._input_arguments['Grading']

        # compute downward force of the foundation
        # only the foundation directly below the caisson
        return (Grading.rho - rho_w)*9.81*(V/H * hb**2 + B*hb)

    def _effective_width_foundation(self, rho_w, hb, B, m):
        """ Compute the effective width of the foundation """
        # slope
        V, H = 1, 1

        # compute moment around the middle of the foundation
        Mb = self.goda.Ma() + self.goda.P() * hb

        # compute the eccentricity
        e = Mb/(self.goda._dFv(m) + self._Fsill(rho_w, hb, B))

        # compute effective width of the foundation
        return B + 2*V/H*hb - 2*e

    def _validate_variant(self, variants):
        """ Validate the input of the variant

        Parameters
        ----------
        variants : tuple
            variantIDs given as args input
        """
        # check if a variant is specified
        if not variants:
            # no input is given so get the valid args for variant
            valid_args = self.variantIDs
            valid_args.append('all')

            # raise error
            valid = ', '.join(valid_args)
            raise InputError(
                'did not specify which variants to use, possible arguments '
               f'are {valid}')

        # check if input all is in variants
        if 'all' in variants:
            # set specified variants to all variantIDs
            variants = tuple(self.variantIDs)

        # return the variants
        return variants

    def _layers(self, variantID):
        """ compute the coordinates of all layers

        Parameters
        ----------
        variantID : str
            identifier of the variant, see :py:attr:`variantIDs` for a
            list of all generated variants.

        Returns
        -------
        dict
            coordinates of all layers
        """
        # get the structure of the current variant
        structure = self.get_variant(variantID)

        # set empty dict to store coordinates in
        coordinates = {}

        # get the slope
        V, H = self._input_arguments['slope_foundation']

        # conmpute constant to transfrom thickness of layer to x and
        # y coordinates, switched V and H because orthogonality
        transform_x = V/np.sqrt(V**2+H**2)
        transform_y = H/np.sqrt(V**2+H**2)

        # get geometrical parameters of the caisson
        B = structure['caisson']['B']
        hc = structure['caisson']['h_acc'] + structure['caisson']['Rc']
        hb = structure['caisson']['hb']
        Bm = structure['caisson']['Bm']

        # determine thickness of the armour layer
        layers_armour = structure['armour']['layers']
        kt_armour = layer_coefficient(
            self._input_arguments['armour'], layers=layers_armour,
            placement='standard')
        dn = structure['armour']['class Dn50']
        t_armour = kt_armour*layers_armour*dn

        # check if scour protection is required
        if Bm >= self.width_scour:
            # no scour protection required
            t_scour = 0
            req_width = 0

        else:
            # define thickness of the scour
            t_scour = 0.4

            # compute required width of the scour protection
            req_width = self.width_scour - Bm - H/V * (hb + t_armour - t_scour)

            if req_width <= 0:
                # negative or zero width required width of the foundation
                # is enough, thus set values to 0
                t_scour, req_width = 0, 0

        # define line of the caisson
        # clockwise starting at the lower left of the caisson
        coordinates['caisson'] = {
            'x': [-0.5*B, -0.5*B, 0.5*B, 0.5*B, -0.5*B],
            'y': [hb, hb + hc, hb + hc, hb, hb]
        }

        # define points of the foundation
        arm_top = hb + t_armour

        arm_x3 = -0.5*B
        arm_x2 = arm_x3 - Bm
        arm_x1 = arm_x2 - H/V * (arm_top - t_scour)

        found_x0 = arm_x1 - req_width - H*t_scour/V
        found_x1 = arm_x1 - req_width
        found_x2 = (arm_x1 + H*(t_armour * transform_y)/V
                       + t_armour*transform_x)
        found_x3 = found_x2 + H/V * (hb - t_scour)
        found_x4 = -0.5*B
        found_x5 = 0.5*B + Bm
        found_x6 = found_x5 + H/V * hb

        # define line of the foundation layers
        # clockwise starting at the left with the upper line
        coordinates['armour'] = {
            'x': [arm_x1, arm_x2, arm_x3, arm_x3, found_x3, found_x2, arm_x1],
            'y': [t_scour, arm_top, arm_top, hb, hb, t_scour, t_scour]
        }

        # check if scour protection is needed
        if req_width != 0:
            # protection needed
            coordinates['foundation'] = {
                'x': [found_x0, found_x1, found_x2, found_x3, found_x4,
                      found_x5, found_x6, found_x0],
                'y': [0, t_scour, t_scour, hb, hb,
                      hb, 0, 0]
            }
        else:
            # not needed
            coordinates['foundation'] = {
                'x': [found_x2, found_x3, found_x4, found_x5, found_x6,
                      found_x2],
                'y': [0, hb, hb, hb, 0,
                      0]
            }

        return coordinates

    def print_logger(self, level='warnings'):
        """ Print messages and warnings in the logger

        Parameters
        ----------
        msg_level : {'info', 'warnings'}, optional, default: 'warnings'
            specify print level, highest level is warnings and lowest
            level is info. Note that the info level will also print all
            warnings
        """
        # check if correct input has been given
        if level.lower() not in ['warnings', 'info']:
            raise NotSupportedError(
                f'{level} not implemented, must be info or warnings')

        # print logger
        for type, messages in self.logger.items():
            if level.lower() == 'warnings':
                if type == 'INFO':
                    continue
            print(f'{type}:')
            if messages:
                for message in messages:
                    print(message)
            else:
                print(f'no {type.lower()} messages in log')
            print()

    def get_variant(self, variantID):
        """ Get the dimensions for the specified variant

        Parameters
        ----------
        variantID : str
            identifier of the variant, see :py:attr:`variantIDs` for a
            list of all generated variants.

        Returns
        -------
        dict
            Parameters and values of the caisson (B, Rc) and the
            foundation layer (Dn50, Rock class) for one variant.

        Raises
        ------
        KeyError
            If there is no variant with the given identifier
        """
        variant = {}
        if variantID in self.variantIDs:
            key_u = self.variantIDs.index(variantID)
        else:
            raise KeyError(f'Variant with ID = {variantID} is not a variant, '
                           f'generated variants are: {self.variantIDs}')

        # add caisson
        variant['caisson'] = self.structure['caisson']

        # add foundation armour and underlayer
        variant['armour'] = self.structure['armour']

        if 'foundation' in self.structure:
            variant['foundation'] = {}
            for param, val in self.structure['foundation'].items():
                if param == 'state':
                    variant['foundation'][param] = val
                else:
                    variant['foundation'][param] = val[key_u]

        # return the generated variant
        return variant

    def print_variant(self, *variants, decimals=3):
        """ Print the details for the specified variant(s)

        This method will print the details of the structure for
        the specified variant(s). It prints the dimensions of the
        caisson (B, Rc) and the dimensions of the foundation layer
        (Dn50, rock class). Furthermore, from the table the normative
        LimitState for the design can be read.

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        decimals : int, optional, default: 3
            number of decimal places to round to

        Raises
        ------
        InputError
            If no arguments are specified
        KeyError
            If there is no variant with the given identifier
        """
        # validate variants
        variants = self._validate_variant(variants)

        # print variant
        for id in variants:
            variant = self.get_variant(id)

            # print the name of the table
            if isinstance(self.id, int):
                table_name = f'Variant {self.id}{id}:'
            else:
                table_name = f'Variant {id}:'

            # set defaults and generate empty tables
            headers_caisson = ['parameter', 'value']
            table_caisson, table_f = [], []
            set_headers = True

            # fill the empty tables
            for i, (layer, dimensions) in enumerate(variant.items()):
                # caisson table
                if layer == 'caisson':
                    for param, val in dimensions.items():
                        if param == 'Rc':
                            val = np.round(val, decimals)
                            state = self.structure['caisson']['state_overtop']
                            label = self._LimitStates[state].label
                            table_val = f'{val} (with {label})'
                            table_caisson.append([param, table_val])
                        elif param == 'B':
                            val = np.round(val, decimals)
                            state = self.structure['caisson']['state_goda']
                            label = self._LimitStates[state].label
                            table_val = f'{val} (with {label})'
                            table_caisson.append([param, table_val])
                        elif param == 'state_overtop' or param == 'state_goda':
                            continue
                        else:
                            table_caisson.append(
                                [param, np.round(val, decimals)])
                # table for the foundation layer
                else:
                    if set_headers:
                        # set headers in first iteration
                        headers_foundation = list(dimensions.keys())
                        headers_foundation.insert(0, 'layer')
                        set_headers = False

                    # make row with the layer and dimensions
                    row = [layer]
                    row.extend(dimensions.values())

                    # check if state is in the headers
                    if 'state' in headers_foundation:
                        # get the normative state and index of state
                        state = dimensions['state']
                        i_state = headers_foundation.index('state')

                        # check if state is an int
                        if isinstance(state, int):
                            # get the label of the normative limitstate
                            label = self._LimitStates[state].label
                            row[i_state] = label

                    # add row to the table
                    table_f.append(row)

            # print tables
            print(table_name)
            print()
            print('  Caisson dimensions')
            print(tabulate(table_caisson, headers_caisson, tablefmt="github"))
            print()
            print('  Foundation dimensions')
            print(tabulate(
                table_f, headers_foundation, tablefmt="github",
                floatfmt=(f'.{decimals}f')))
            print('\n')

    def plot(self, *variants, wlev=None, save_name=None):
        """ Plot the cross section of the specified breakwater(s)

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        wlev : str, optional, default: None
            label of the :py:class:`LimitState` from which the water
            level will be plotted. If no value is specified the water
            level from the normative limit state is used, which is the
            normative LimitState from the crest freeboard computation.
        save_name : str, optional, default: None
            if given the cross section is not shown but saved with the
            given name

        Raises
        ------
        InputError
            If no variants are specified or if the label of wlev is not
            a valid label of a :py:class:`LimitState`
        KeyError
            If there is no variant with the given identifier
        """
        # validate variants
        variants = self._validate_variant(variants)

        if wlev is None:
            wlev = self.structure['caisson']['state_overtop']
        else:
            for i, LimitState in enumerate(self._LimitStates):
                if LimitState.label == wlev:
                    wlev = i
                    break

        # check if wlev has been changed from LimitState label into index
        if isinstance(wlev, str):
            # wlev is still a string so not changed, which means that
            # the specified wlev is not a specified LimitState
            raise InputError('There is no LimitState with the given label')

        # set figure
        plt.figure(figsize=(10, 5))

        for i, id in enumerate(variants):
            if len(variants) == 2:
                # make subplot if two variants must be plotted
                plt.subplot(1, 2, i+1)

            # get the coordinates
            coordinates = self._layers(id)

            # set xlim_max and xlim_min variable
            xlim_max, xlim_min = 0, 0

            # plot lines
            for layer, lines in coordinates.items():
                plt.plot(lines['x'], lines['y'], color='k')

                # check largest value for xlim
                if np.max(lines['x']) >= xlim_max:
                    # set max as xlim_max
                    xlim_max = np.max(lines['x'])

                # check smallest value for xlim
                if np.min(lines['x']) <= xlim_min:
                    # set min as xlim_min
                    xlim_min = np.min(lines['x'])

            # plot bottom and wlev
            x_wlev_max = np.max(coordinates['caisson']['x'])

            plt.axhline(y=0, color='k', linewidth=2)
            plt.hlines(
                y=self._LimitStates[wlev].h, xmin=xlim_min*1.2,
                xmax=-x_wlev_max, color='b')
            plt.hlines(
                y=self._LimitStates[wlev].h, xmin=x_wlev_max,
                xmax=xlim_max*1.2, color='b')

            # set xlim and ylim
            ymax = np.max(coordinates['caisson']['y'])*1.2
            plt.xlim(xlim_min*1.2, xlim_max*1.2)
            plt.ylim(-0.5, ymax)

            # add title to the plot
            if save_name is None:
                if isinstance(self.id, int):
                    title = ('Cross section of monolithic breakwater '
                             f'{self.id}{id}')
                else:
                    title = f'Cross section of monolithic breakwater {id}'
            else:
                name = save_name.split('/')[-1]
                title = f'Cross section of {name}'

            # add title, grid and set equal aspect ratio
            plt.title(title)
            plt.grid()
            plt.gca().set_aspect('equal', adjustable='box')

        plt.tight_layout()

        # save the figure
        if save_name is not None:
            plt.savefig(f'{save_name}.png')
            plt.close()
        else:
            plt.show()

    def plot_pressure(self):
        """ Plot pressure distribution computed extended Goda formula

        Plots the pressure distribution computed with the extended
        Goda formula (Takahasi, 2002) together with the dimensions
        of the monolithic breakwater.

        .. warning::
           Do not read the dimensions of the monolithic breakwater from
           the axes of the figure. The correct dimensions of the
           monolithic breakwater are depicted in the figure, or use
           :py:meth:`get_variant` or :py:meth:`print_variant`.
        """
        self.goda.plot()

    def area(self, variantID):
        """ Compute the area of all layers

        Method computes the area of each layer using Gauss's area
        formula. Which is given by the following formula:

        .. math::
           \\mathbf{A}=\\frac{1}{2} | \\sum_{i=1}^{n-1} x_{i} y_{i+1}
           +x_{n} y_{1}-\\sum_{i=1}^{n-1} x_{i+1} y_{i}-x_{1} y_{n} |

        Parameters
        ----------
        variantID : str
            identifier of the variant, see :py:attr:`variantIDs` for a
            list of all generated variants.

        Returns
        -------
        dict
            dict with the area of each layer
        """
        # get the coordinates of the layers
        coordinates = self._layers(variantID)

        # iterate over the layers to compute the area
        area = {}
        for layer, coord in coordinates.items():
            # get the x and y coordinates
            x = coord['x']
            y = coord['y']

            # use Gauss's area formula
            A = 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

            # add to area dict
            area[layer] = A

        return area

    def cost(
            self, *variants, concrete_price, fill_price, unit_price=None,
            output='variant'):
        """ Compute the cost per meter for each variant

        Method to compute the cost of each generated variant, the cost
        is computed per meter

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        concrete_price : float
            price of concrete per m³
        fill_price : float
            price of the fill material per m³
        unit_price : float, optional, default: None
            the cost of an armour unit per m³, required if the armour
            layer of the foundation is made out of armour units
        output : {variant, layer, average}
            format of the output dict, variant returns the total cost
            of each variant, layer the cost of each layer for each
            variant and average returns the average cost.

        Returns
        -------
        dict
            the cost

        Raises
        ------
        RockGradingError
            if no pricing is included in the given RockGrading
        """
        # validate variants
        variants = self._validate_variant(variants)

        # get the grading, and check if the cost has been added
        Grading = self._input_arguments['Grading']

        if 'price' in Grading[list(Grading.grading.keys())[0]]:
            # pricing has been added
            pass
        else:
            # pricing has not been added, raise error
            raise RockGradingError('There is no pricing in the RockGrading')

        # set empty dict to store the output in
        cost = {}

        # iterate over the generated variants
        for id in variants:
            # get the areas and structure of the variants
            areas = self.area(id)
            structure = self.get_variant(id)

            # iterate over the layers to price each layer
            variant_price = {}
            for layer, area in areas.items():
                if layer is 'caisson':
                    # compute price of the caisson
                    Pc = structure['caisson']['Pc']
                    price = area*Pc*concrete_price + area*(1-Pc)*fill_price
                elif (self._input_arguments['armour'] is not 'Rock'
                        and layer is 'armour'):
                    # concrete armour units
                    if unit_price is not None:
                        price = area * unit_price
                    else:
                        raise InputError(
                            ('argument unit_price is required when computing '
                             'the cost with armour units as BermMaterial'))
                else:
                    # layer of the breakwater
                    rock_class = structure[layer]['class']

                    # get the price per meter
                    price = Grading[rock_class]['price'] * area

                # add to dict
                variant_price[layer] = np.round(price, 2)

            # add to cost dict
            if output is 'variant' or output is 'average':
                # add total cost of all layers
                cost[id] = np.round(np.sum(list(variant_price.values())), 2)
            elif output is 'layer':
                # add the cost of each layer
                cost[id] = variant_price
            else:
                # invalid input
                raise NotSupportedError(
                    (f'Cost can\'t be exported as {output}, must be variant, '
                      'layer or average'))

        # check if average must be computed
        if output is 'average':
            # compute average cost
            cost = {'average': np.round(np.average(list(cost.values())), 2)}

        # check if dry dock has been added
        if self.price is None:
            # no dry dock added
            self.price = cost
        else:
            # check output
            if output == 'variant':
                # add investment to each variant
                for id, price in cost.items():
                    cost[id] = price + self.price

            elif output == 'layer':
                # add investment as a key
                cost['investment'] = self.price

            elif output == 'average':
                # add investment
                cost['average'] += self.price

        # check if the average cost is computed for only one variant
        if len(variants) == 1 and output == 'average':
            # print user_warning and change key into variant
            cost[variants[0]] = cost.pop('average')
            user_warning(
                ('Computing the average for one variantID, changed key '
                 'average in dict with the specified variantID'))

        return cost

    def dry_dock(self, investment, length):
        """ Add the investment cost of a dry dock to the concept

        This method adds the investment required to rent a dry dock to
        the concept. The investment cost is added to the concept by
        dividing the investment trough the length of the breakwater to
        get the required investment per running meter.

        Parameters
        ----------
        investment : float
            the investment required to rent a dry dock
        length : float
            length of the breakwater [m]
        """
        # check if cost have already been added
        if self.price is not None:
            # iterate over the price dict
            for id, price in self.price.items():
                # check if mode was average
                if id == 'average':
                    # mode was average
                    self.price[id] = price + investment/length
                    break
                else:
                    # check if price is a dict
                    if isinstance(price, dict):
                        # mode was layer
                        # add investment as additional key
                        self.price['investment'] = investment/length
                        break
                    else:
                        # set new price
                        self.price[id] = price + investment/length

        else:
            # set investment as attribute
            self.price = investment/length
