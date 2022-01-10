import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import pandas as pd
import os
from warnings import catch_warnings
from tabulate import  tabulate
import pickle

from .utils._kwarg_validator import _process_kwargs, _RM_vkwargs, _C_vkwargs
from .utils._design_explorer import _DE_params
from .utils.exceptions import RockGradingError, ArmourUnitsError, InputError, user_warning, NotSupportedError
from .utils._progress import ProgressBar
from .utils.cost import _process_cost
from .conditions import LimitState


from breakwater.shape_3D.shape import all_designs, wave_angles_structure
from breakwater.utils.exceptions import Design3DError

class structure_3D():

    def __init__(
                self,
                kml_path,
                wave_conditions,
                slope,
                slope_foreshore,
                B,
                N,
                rho_w,
                ArmourUnit,
                Grading,
                core_material,
                LimitState,
                shape = 'LinearRing',
                structure_type= 'breakwater',
                wave_direction= 'right',
                safety=1,
                slope_toe=(2, 3),
                B_toe=None,
                layers=1,
                layers_underlayer=2,
                vdm= 'max',
                filter_rule=None,
                Soil=None,
                phi=40,
                id=None,
                **kwargs
                ):
        """
        Class structure_3D which can be used to expand the 2D design to a 3D structure

        Parameters
        ----------
        kml_path: str
            path with the kml file
        shape: str
            shape of the structure. default is LinearRing
        wave_conditions: dict
            dictionary with orientation the key and item is a dict with Hm0 and Tp
        structure_type: str, default is breakwater
            Either revetment or breakwater
        wave_direction: str, default= 'right'
            When walking in clockwise direction along the coordinates, do the waves
            come from the right or the left
        h: float
            Water depth at the structure
        Sd: int
            Allowed damage level of the armour
        Nod: int
            Allowed damage level at the toe
        q: int
            Allowed number of overtopping in q/l/m
        slope : tuple
            Slope of the armour layer (V, H). For example a slope of 3V:4H
            is defined as (3, 4)
        slope_foreshore : tuple
            slope of the foreshore (V, H). For example a slope of 1:100 is
            defined as (1, 100)
        rho_w : float
            density of water [kg/m³]
        B : float
            Crest width [m]
        N : int
            Number of incident waves at the toe of the structure [-]
        LimitState : :py:class:`LimitState` or list of :py:class:`LimitState`
            ULS, SLS or another limit state defined with
            :py:class:`LimitState`
        Grading : :py:class:`RockGrading`
            standard rock grading defined in the NEN-EN 13383-1 or a user
            defined rock grading
        core_material : dict
            core_class and Dn50
        safety : float, optional, default: 1
            safety factor of design (number of standard deviations from the
            mean)
        slope_toe : tuple, optional, default: (2,3)
            slope of the toe
        B_toe : float, optional, default: None
            width of the top of the toe in meters. By default the width of
            toe is taken as 3 * Dn50_toe.
        beta : float, optional, default: 0
            angle between direction of wave approach and a line normal to
            the breakwater (degrees).
        layers : int, optional, default: 2
            number of layers in the armour layer
        layers_underlayer : int, optional, default: 2
            number of layers in the underlayer
        vdm : {min, max, avg}, optional, default: max
            value to return in case both the deep and shallow water formula
            are valid. min for the lowest value, max for the highest value
            and avg for the average value, default is max.
        Soil : :py:class:`Soil`, optional, default: None
            by default Soil is None, which means that the geotechnical checks
            are not performed. By specifying a Soil object, the geotechnical
            checks are automatically performed.
        phi : float, optional, default: 40
            internal friction angle of rock [degrees]
        id : int, optional, default: None
            add a unique id to the breakwater
        """

        self.kml_path = kml_path
        self.wave_conditions = wave_conditions
        self.shape = shape
        self.structure_type = structure_type
        self.wave_direction = wave_direction
        self.Grading = Grading

        LimitState = LimitState.Limit_states

        designs_all_directions = all_designs(
                                            kml_path= self.kml_path,
                                            wave_direction= self.wave_direction,
                                            wave_conditions= self.wave_conditions,
                                            shape= self.shape,
                                            structure_type= self.structure_type,
                                            slope= slope,
                                            slope_foreshore= slope_foreshore,
                                            B= B,
                                            N= N,
                                            rho_w= rho_w,
                                            LimitState= LimitState,
                                            ArmourUnit= ArmourUnit,
                                            Grading= Grading,
                                            core_material= core_material,
                                            safety= safety,
                                            slope_toe= slope_toe,
                                            B_toe= B_toe,
                                            layers= layers,
                                            layers_underlayer= layers_underlayer,
                                            vdm= vdm,
                                            filter_rule= filter_rule,
                                            Soil= Soil,
                                            phi= phi,
                                            id= id,
                                            **kwargs)

        for part, specs in designs_all_directions.items():
            for i in range(len(specs['wave_conditions'])):
                condition = specs['wave_conditions'][i]
                c = 0
                for design in condition['design']:
                    if len(designs_all_directions[part]['coordinates'][c]) == 2:
                        designs_all_directions[part]['coordinates'][c].append(float('-inf'))
                    if design != None:
                        if type(designs_all_directions[part]['coordinates'][c][-1]) == float:
                            if design.structure['armour']['class Dn50'] > designs_all_directions[part]['coordinates'][c][-1] and \
                               design.Rc > designs_all_directions[part]['coordinates'][c][-1]:
                                designs_all_directions[part]['coordinates'][c][-1] = design
                        else:
                            if design.structure['armour']['class Dn50'] > \
                                    designs_all_directions[part]['coordinates'][c][-1].structure['armour']['computed Dn50'] and \
                                    design.Rc > designs_all_directions[part]['coordinates'][c][-1].Rc:
                                designs_all_directions[part]['coordinates'][c][-1] = design
                    c += 1

        d = {
        (i, j): [designs_all_directions[i]['coordinates'][j][:-1],designs_all_directions[i]["coordinates"][j][-1], designs_all_directions[i]["distance"][j]]
        for i in designs_all_directions.keys()
        for j in range(len(designs_all_directions[i]['coordinates']))
            }

        mux = pd.MultiIndex.from_tuples(d.keys())
        df = pd.DataFrame(list(d.values()), index=mux)

        df.rename(columns= {0: 'coordinates', 1: 'structure', 2: 'length [m]'}, inplace= True)
        df['class Dn50'] = [None] * len(df)
        df['Rc'] = [None] * len(df)


        i = 0
        for index, row in df.iterrows():
            df['class Dn50'].iloc[i] = row['structure'].structure['armour']['class Dn50']
            df['Rc'].iloc[i] = row['structure'].Rc
            i += 1

        sections = df.index.get_level_values(0).unique()
        new_df = pd.DataFrame()

        for s in sections:
            df_sec = df.xs(s).sort_values(by = ['class Dn50', 'Rc'], ascending= False).reset_index(drop= True)
            length = sum(df_sec['length [m]'])
            df_sec['length [m]'] = [length] * len(df_sec)

            first_row = df_sec.loc[0]
            if len(df_sec) > 1:
                for i in range(0, len(df.xs(s))):
                    other_row = df.xs(s).loc[i]
                    if i == 0:
                        first_row['coordinates'] = other_row['coordinates']
                    else:
                        first_row['coordinates'].extend(other_row['coordinates'])
            new_df = new_df.append(first_row)

        new_df['sections'] = sections
        new_df = new_df.set_index('sections', drop= True)

        self.df_design = new_df

    def plot_topview(self):
        """
        Plot the top view of the structure shape

        Returns
        -------
        Matplotlib.figure

        """

        orientation_dict = wave_angles_structure(kml_path= self.kml_path, wave_conditions= self.wave_conditions,
                                                   wave_direction = self.wave_direction, shape = self.shape)

        color = iter(cm.tab20b(np.linspace(0, 1, len(orientation_dict))))
        plt.figure(figsize= (10, 6.7))
        for part, specs in orientation_dict.items():
            label = True
            colr = next(color)
            for c in specs['coordinates']:
                if label == True:
                    plt.plot([c[0][0], c[1][0]], [c[0][1], c[1][1]], label= part, color = colr)
                    label = False
                else:
                    plt.plot([c[0][0], c[1][0]], [c[0][1], c[1][1]], color= colr)

        plt.legend()
        plt.axis('off')


    def plot_section(
                    self,
                    *variants,
                    section,
                    wlev=None,
                    save_name=None,
                    ):
        """Plot the cross section of the specified section

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        section: tuple
            section to be plotted
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

        if section not in self.df_design.index.tolist():
            raise Design3DError(f'{section} does not exist. Choose one of {self.df_design.index.get_level_values(0).unique().values}')

        section = self.df_design.loc[section]

        section['structure'].plot(*variants, wlev= wlev, save_name= save_name)

    def print_section(self, *variants, section, decimals=3):
        """Print the details for the specified variant(s)

        This method will print the computed Dn50, rock class, average
        Dn50 of the class, normative LimitState for all layers and
        specified variants.

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        section: tuple
            Section to be printed
        decimals : int, optional, default: 3
            number of decimal places to round to

        Raises
        ------
        InputError
            If no arguments are specified
        KeyError
            If there is no variant with the given identifier
        """
        if section not in self.df_design.index.tolist():
            raise Design3DError(f'{section} does not exist. Choose one of {self.df_design.index.get_level_values(0).unique().values}')

        section = self.df_design.loc[section]

        section['structure'].print_variant(*variants, decimals= decimals)

    def plot_3D(self):

        plt.figure(figsize= (10, 6.7))

        start = None

        for index, row in self.df_design.iterrows():
            if index == 'section A':
                coords = row['coordinates']
                lon, lat = list(zip(*coords))
                if start == None:
                    start = (lon[0], lat[0])

                lon = np.array(lon) - start[0]
                lat = np.array(lat) - start[1]

                coords = row['structure'].plot('a', get_data = True)
                for layer, line in coords.items():
                    x = line['x']
                    y = line['y']
                    if self.wave_direction == 'right':
                        x = np.array(x) * -1
                    plt.plot(x, y, color= 'b')
                    plt.gca().set_aspect("equal", adjustable="box")


    def totalcost_3D(
                    self,
                    core_price,
                    grading_price,
                    unit_price=None,
                    equipment=None,
                    transport_cost= None,
                    algorithm = 'smart_combinations',
                    limit = None,
                    optimize_on = ['cost', 'time'],
                    output="variant",
                    ):
        """

        Parameters
        ----------
        core_price : dict
            cost of the core material per m³
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        unit_price : dict
            the cost of an armour unit per m³
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        grading_price: dict
            cost of material and CO2 per stone grading. {'LMA_40/60': {'cost': ..., 'CO2':...}}
        equipment: lst
            list of equipment out of Equipment class
        transport_cost : dict
            the cost to transport a m³ of rock from the quarry to the
            project location
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        algorithm: str
            Which algorithm should be used to determine the set of equipment
            either smart_combinations, cheap_combinations. Default is smart_combinations
        limit: int
            In the case the cheap_combinations algorithm is chosen this limit gives the maximum number
            of equipment used
        output : {variant, layer, average}
            format of the output dict, variant returns the total cost
            of each variant, layer the cost of each layer for each
            variant and average returns the average cost.

        Returns
        -------
        tuple
            the cost

        Raises
        ------
        RockGradingError
            if no pricing is included in the given RockGrading
        EquipmentError
            if the given equipment is not able to fill all the sections

        """

        cost_df = self.df_design.copy()

        self.Grading.add_cost(grading_price)

        lst_total_cost = []
        lst_total_CO2 = []
        lst_duration = []
        lst_opt_equip = []

        for index, row in self.df_design.iterrows():
            variants = row['structure']
            length = row['length [m]']

            total_cost, total_CO2, duration, opt_equip = variants.total_cost(
                                                                                'all',
                                                                                core_price= core_price,
                                                                                unit_price= unit_price,
                                                                                equipment= equipment,
                                                                                length= length,
                                                                                plot_error= False,
                                                                                transport_cost= transport_cost,
                                                                                algorithm = algorithm,
                                                                                limit = limit,
                                                                                optimize_on = optimize_on,
                                                                                output= output,
                                                                                )

            lst_total_cost.append(total_cost)
            lst_total_CO2.append(total_CO2)
            lst_duration.append(duration)
            lst_opt_equip.append(opt_equip)

        if equipment != None:
            cost_df['total cost [EUR]'] = lst_total_cost
            cost_df['total CO2 [kge]'] = lst_total_CO2
            cost_df['duration [wk]'] = lst_duration
            cost_df['equipment'] = lst_opt_equip
        else:
            cost_df['material cost [EUR]'] = lst_total_cost
            cost_df['material CO2 [kge]'] = lst_total_CO2

        self.df = cost_df

        return cost_df

class Configurations_3D:

    def __init__(self,
                 kml_path,
                 wave_conditions,
                 structure,
                 LimitState,
                 rho_w,
                 slope_foreshore,
                 Grading,
                 N,
                 shape = 'LinearRing',
                 wave_direction= 'right',
                 Soil=None,
                 safety=1,
                 **kwargs):

        data = wave_angles_structure(kml_path, wave_conditions, LimitStates= LimitState.Limit_states, wave_direction = 'right', shape = 'Linestring')

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

        # if material cost is added the standard names are below
        # if also equipment is added it is changed to total_cost and total_CO2 (see add_cost)
        self.col_cost = 'material_cost'
        self.col_CO2 = 'material_CO2'


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
            core_material = concept['core_material']
            B_toe = concept['B_toe']
            slope = concept['slope']
            slope_toe = concept['slope_toe']

            if RRM_compute:

                try:
                    with catch_warnings(record=True) as w:
                        RM_rock = structure_3D(
                            kml_path= kml_path, wave_direction= wave_direction,
                            wave_conditions= wave_conditions, shape= shape,
                            structure_type= 'breakwater', slope=slope, slope_foreshore=slope_foreshore,
                            rho_w=rho_w, B=B, N=N, LimitState=LimitState,
                            Grading=Grading, core_material= core_material,
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
                                             'Dn50_core': [core_material['Dn50']],
                                             'B_toe': [B_toe],
                                             'slope': [slope],
                                             'slope_toe': [slope_toe],
                                             'warnings': [w]})
                self.df = self.df.append(temp_df, ignore_index=True, sort=True)

                RM_bar.next()

            if CRM_compute:
                try:
                    with catch_warnings(record=True) as w:
                        RM_units = structure_3D(
                            kml_path= kml_path, wave_direction= wave_direction,
                            wave_conditions= wave_conditions, shape= shape, N= N,
                            structure_type= 'breakwater', slope=slope,
                            slope_foreshore=slope_foreshore, B=B,
                            rho_w=rho_w, LimitState=LimitState, safety=safety,
                            Grading=Grading, ArmourUnit=ArmourUnit, phi=phi,
                            core_material= core_material, slope_toe=slope_toe,
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
                                             'Dn50_core': [core_material['Dn50']],
                                             'B_toe': [B_toe],
                                             'slope': [slope],
                                             'slope_toe': [slope_toe],
                                             'warnings': [w]})
                self.df = self.df.append(temp_df, ignore_index=True, sort=True)

                RM_bar.next()

            if CRMR_compute:
                try:
                    with catch_warnings(record=True) as w:
                        RM_units = structure_3D(
                            kml_path= kml_path, wave_direction= wave_direction,
                            wave_conditions= wave_conditions, shape= shape,
                            structure_type= 'revetment', slope=slope,
                            slope_foreshore=slope_foreshore, B=B, N= N,
                            rho_w=rho_w, LimitState=LimitState, safety=safety,
                            Grading=Grading, ArmourUnit=ArmourUnit, phi=phi,
                            core_material= core_material, slope_toe=slope_toe,
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
                                             'Dn50_core': [core_material['Dn50']],
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
            self,
            core_price,
            grading_price,
            unit_price = None,
            equipment=None,
            transport_cost= None,
            algorithm = 'smart_combinations',
            limit = None,
            optimize_on = ['cost', 'time'],
            concrete_price= None,
            fill_price= None,
            output="variant"):
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
        optimize_on: str/list
            On what to optimize
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
            _process_cost(structure, cost, self._Grading)

        # set list to store cost in

        # iterate over the generated concepts
        for i, row in self.df.iterrows():
            # check if concept is not None
            if row.concept is None:
                # not a valid concept
                continue
            else:
                # valid concept
                # check types and compute price

                if row.type == 'RRM':
                    df_cost = row.concept.totalcost_3D(
                        core_price= core_price, unit_price= unit_price, grading_price= grading_price,
                        equipment = equipment, transport_cost=transport_cost, algorithm= algorithm,
                        limit= limit, optimize_on= optimize_on, output= output)


                elif row.type == 'CRM':
                    df_cost  = row.concept.totalcost_3D(
                        core_price= core_price, unit_price= unit_price, grading_price= grading_price,
                        equipment = equipment, transport_cost=transport_cost, algorithm= algorithm,
                        limit= limit, optimize_on= optimize_on, output= output)

                elif row.type == 'CRMR':
                    df_cost = row.concept.totalcost_3D(
                        core_price= core_price, unit_price= unit_price, grading_price= grading_price,
                        equipment = equipment, transport_cost=transport_cost, algorithm= algorithm,
                        limit= limit, optimize_on= optimize_on, output= output)
                    print(row.concept.df_design)
                elif row.type == 'RC' or row.type == 'CC':
                    # check if investment cost must be added

                    df_cost  = row.concept.totalcost_3D(
                        core_price= core_price, unit_price= unit_price, grading_price= grading_price,
                        equipment = equipment, transport_cost=transport_cost, algorithm= algorithm,
                        limit= limit, optimize_on= optimize_on, output= output)

                else:
                    raise NotSupportedError(f'{row.type} is not supported')

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
        | cost                     |     o      |     o      |
        +----------------------------------------------------+
        | CO2                      |     o      |     o      |
        +----------------------------------------------------+
        | install_duration         |     o      |     o      |
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

        if 'cost' in params:
           params[params.index('cost')] = self.col_cost
        if 'CO2' in params:
           params[params.index('CO2')] = self.col_CO2

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
                if row[self.col_cost][id] is not None:
                    file_name = f'{row.type}.{row.id}'
                    if num_concepts > 1:
                        # add id if more than 1 concept
                        file_name = f'{file_name}{id}'

                    save_name = f'{new_dir}/{file_name}'

                    # save cross section of the current concept
                    # add equipment to plot if present
                    if 'optimal_equipment' in row.index.values:
                        row.concept.plot(id, save_name=save_name, equipment= row.optimal_equipment[id])
                        # add equipment to data
                    else:
                        row.concept.plot(id, save_name=save_name)

                    data = {'CaseNo': [CaseNo], 'type': [row.type]}

                    # get the values from the variant and add to data
                    variant = row.concept.get_variant(variantID=id)
                    to_explorer = _DE_params(
                        args=params, variant=variant, row=row, concept=row.concept,
                        structure=row.type, slopes=slopes,
                        change_CRM_class=format_class_as_string)
                    data.update(to_explorer)

                    # check if cost must be included
                    if self.col_cost in params:
                        # check if cost have been added
                        if self.col_cost in self.df.columns:
                            # add cost to data
                            data[self.col_cost] = round(row[self.col_cost][id], 2)

                        else:
                            raise KeyError(
                                'cost have not been added, use add_cost to add '
                                'EUR cost to the df')

                    if self.col_CO2 in params:
                        if self.col_CO2 in self.df.columns:
                            # add cost to data
                            data[self.col_CO2] = round(row[self.col_CO2][id], 4)

                        else:
                            raise KeyError(
                                'CO2 have not been added, use add_cost to add CO2'
                                'cost to the df')

                    if 'install_duration' in params:
                        if 'install_duration' in self.df.columns:
                            # add cost to data
                            data['install_duration'] = round(row.install_duration[id], 4)

                        else:
                            raise KeyError(
                                'equipment have not been added, thus no installation rates are provided')

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

