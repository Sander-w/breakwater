import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import pandas as pd

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
            cost_df['total cost'] = lst_total_cost
            cost_df['total CO2'] = lst_total_CO2
            cost_df['duration'] = lst_duration
            cost_df['equipment'] = lst_opt_equip
        else:
            cost_df['material_cost'] = lst_total_cost
            cost_df['material_CO2'] = lst_total_CO2

        self.df = cost_df

        return cost_df

