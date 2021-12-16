import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
import pandas as pd

from breakwater.shape_3D.shape import all_designs, wave_angles_structure
from breakwater.utils.exceptions import Design3DError


class structure_3D():

    def __init__(
                self,
                kml_path,
                wave_conditions,
                shape = 'LinearRing',
                structure_type= 'breakwater',
                wave_direction= 'right',
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
        """

        self.kml_path = kml_path
        self.wave_conditions = wave_conditions
        self.shape = shape
        self.structure_type = structure_type
        self.wave_direction = wave_direction
        self.Grading = None
        self.df_design = pd.DataFrame()

    def plot_topview(self):
        """
        Plot the top view of the structure shape

        Returns
        -------
        Matplotlib.figure

        """

        orientation_dict = wave_angles_structure(kml_path= self.kml_path, wave_conditions= self.wave_conditions,
                                                   wave_direction = self.wave_direction, shape = self.shape)

        color = iter(cm.tab20c(np.linspace(0, 1, len(orientation_dict))))

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

    def design_3D(
                    self,
                    h,
                    Sd,
                    Nod,
                    q,
                    slope,
                    slope_foreshore,
                    B,
                    N,
                    rho_w,
                    ArmourUnit,
                    Grading,
                    core_material,
                    limitstate_label= 'ULS',
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
                    **kwargs):

        """
        Creates a design for each part of the structure by the governing wave direction

        Parameters
        ----------
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

        Returns
        -------
        dict
        """

        self.Grading = Grading

        designs_all_directions = all_designs(
                                            kml_path= self.kml_path,
                                            wave_direction= self.wave_direction,
                                            wave_conditions= self.wave_conditions,
                                            shape= self.shape,
                                            structure_type= self.structure_type,
                                            h= h,
                                            Sd= Sd,
                                            Nod= Nod,
                                            q= q,
                                            slope= slope,
                                            slope_foreshore= slope_foreshore,
                                            B= B,
                                            N= N,
                                            rho_w= rho_w,
                                            ArmourUnit= ArmourUnit,
                                            Grading= Grading,
                                            core_material= core_material,
                                            limitstate_label= limitstate_label,
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
                            if design.structure['armour']['computed Dn50'] > designs_all_directions[part]['coordinates'][c][-1]:
                                designs_all_directions[part]['coordinates'][c][-1] = design
                        else:
                            if design.structure['armour']['computed Dn50'] > designs_all_directions[part]['coordinates'][c][-1].structure['armour']['computed Dn50'] :
                                designs_all_directions[part]['coordinates'][c][-1] = design
                    c += 1

        d = {
        (i, j): [designs_all_directions[i]["coordinates"][j][-1], designs_all_directions[i]["distance"][j]]
        for i in designs_all_directions.keys()
        for j in range(len(designs_all_directions[i]['coordinates']))
            }

        mux = pd.MultiIndex.from_tuples(d.keys())
        df = pd.DataFrame(list(d.values()), index=mux)

        df.rename(columns= {0: 'structure', 1: 'length [m]'}, inplace= True)

        self.df_design = df

        return df

    def plot_section(
                    self,
                    section
                    ):
        """

        Parameters
        ----------
        section: str


        Raises
        ---------
        Design3DError
            If the section does not exist
        """

        if section not in self.df_design.index.tolist():
            raise Design3DError(f'{section} does not exist')

        section = self.df_design.loc[section]

        section['structure'].plot('all')

    def totalcost_3D(
                    self,
                    core_price,
                    unit_price,
                    grading_price,
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
        core_price
        unit_price
        grading_price
        equipment
        transport_cost
        algorithm
        limit
        optimize_on
        output

        Returns
        -------

        """

        if self.df_design.empty:
            raise Design3DError('First create designs for all the different sections of the structure')

        else:
            self.Grading.add_cost(grading_price)

            for index, row in self.df_design.iterrows():
                variants = row['structure']
                length = row['length [m]']

                total_cost, total_CO2, duration, opt_equip = variants.total_cost(
                                                                                    'all',
                                                                                    core_price= core_price,
                                                                                    unit_price= unit_price,
                                                                                    equipment= equipment,
                                                                                    plot_error= False,
                                                                                    transport_cost= transport_cost,
                                                                                    algorithm = algorithm,
                                                                                    limit = limit,
                                                                                    optimize_on = optimize_on,
                                                                                    output= output,
                                                                                    )

                print(total_cost)
                if optimize_on[0] == 'cost':
                    variant = min(total_cost, key = total_cost.get)
                if optimize_on[0] == 'CO2':
                    variant = min(total_CO2, key = total_CO2.get)
                else:
                    variant = min(duration, key = duration)

                total_cost = total_cost[variant] * length
                material_CO2 = material_CO2[variant] * length
                duration = duration[variant] * length
