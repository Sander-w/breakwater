import pandas as pd
import numpy as np
from breakwater.utils.exceptions import InputError

def cheap_combinations(dataframe, optimize_on = 'cost'):

    """
    Parameters
    ----------
    dataframe
    optimize_on: str or list
        one or more of 'CO2', 'cost' or 'time'

    Returns
    -------
    DataFrame
    """

    equipment = dataframe.columns
    equipment = [e for e in equipment if e != 'area']

    if type(optimize_on) == str:
        optimize_on1 = optimize_on
    if type(optimize_on) == list:
        optimize_on1 = optimize_on[0]
        optimize_on2 = optimize_on[1]

    equip_str = ''
    equip_lst = []
    equip_cost = {'cost': 0, 'CO2': 0, 'time': 0}

    for index, row in dataframe.iterrows():
        opt1 = float('inf')
        opt2 = float('inf')

        equip_cheap = None
        dict_optimizations = {}

        for equip in equipment:
            if row[equip] != None:
                if row[equip][optimize_on1] < opt1:
                    opt1 = row[equip][optimize_on1]
                    opt2 = row[equip][optimize_on2]
                    equip_cheap = equip
                    dict_optimizations = row[equip]

                if row[equip][optimize_on1] == opt1 and row[equip][optimize_on2] < opt2:
                    opt1 = row[equip][optimize_on1]
                    opt2 = row[equip][optimize_on2]
                    equip_cheap = equip
                    dict_optimizations = row[equip]


        if equip_cheap not in equip_lst:
            equip_str += (equip_cheap.name + ' + ')
            equip_lst.append(equip_cheap)

        for key, item in dict_optimizations.items():
            equip_cost[key] += item

    equip_combi = {}
    equip_combi[equip_str[:-3]] = equip_cost

    df_optimal = pd.DataFrame.from_dict(equip_combi, orient='index')

    df_optimal['cost'] = df_optimal['cost'].round(2)
    df_optimal['CO2'] = df_optimal['CO2'].round(5)
    df_optimal['time'] = df_optimal['time'].round(5)

    return df_optimal

def cheap_combinations_limit(dataframe, limit, optimize_on = 'cost'):
    """

    Parameters
    ----------
    dataframe
    limit
    optimize_on

    Returns
    -------

    """

    if type(optimize_on) == str:
        optimize_on1 = optimize_on
    if type(optimize_on) == list:
        optimize_on1 = optimize_on[0]
        optimize_on2 = optimize_on[1]

    df_optimal = cheap_combinations(dataframe= dataframe, optimize_on = optimize_on)

    equipment_names = df_optimal.index.tolist()[0].split(' + ')

    for c in [c for c in dataframe.columns if c != 'area']:
        if c.name not in equipment_names:
            del dataframe[c]

    del dataframe['area']

    while len(dataframe.columns) > limit:

        remaining_equipment = [c for c in dataframe.columns]
        marginal_costs = float('inf')
        remove_equip = None
        new_df_optimal = None
        new_dataframe = None

        for equip in remaining_equipment:
            dataframe2 = dataframe.copy()
            del dataframe2[equip]
            dataframe2.dropna(inplace= True, how= 'all')
            if len(dataframe2) == len(dataframe):
               df_optimal2 = cheap_combinations(dataframe2, optimize_on = optimize_on)
               marginal = df_optimal2[optimize_on1].values[0] - df_optimal[optimize_on1].values[0]
               if marginal < marginal_costs:
                   marginal_costs = marginal
                   remove_equip = equip
                   new_df_optimal = df_optimal2
                   new_dataframe = dataframe2

        if remove_equip != None:
            del dataframe[remove_equip]
            df_optimal = new_df_optimal
            dataframe = new_dataframe
        else:
            raise InputError(f'It is not possible to install the structure with only {limit} installation methods. The maximum allowed'
                             f' limit is {len(remaining_equipment)}.')

    return df_optimal


def smart_combinations(dataframe, optimize_on ="cost"):

    """

    Parameters
    ----------
    dataframe
    optimize_on: str or list
        one or more of 'CO2', 'cost' or 'time'

    Returns
    -------
    DataFrame
    """

    df = dataframe

    df.index = np.arange(0, len(df))

    # All the sections
    all_sections = df.index
    # All the equipment
    equipment = [c for c in df.columns if c != 'area']
    equipment_layers = {}
    # Which equipment can install which sections
    for i in range(len(equipment)):
        layers_equip_i = list(df[df[equipment[i]].notna()].index.values)
        equipment_layers[(equipment[i].name, equipment[i].type)] = layers_equip_i


    combinations_dict = []
    combinations_lst = []

    all_equips_dict = {}

    for equip in equipment:

        equip = (equip.name, equip.type)
        i = 0
        installed = equipment_layers[equip].copy()
        equips_dict = {}
        equips_dict[equip[0]] = installed.copy()
        equips_lst = [equip[0]]

        # can't install a single layer
        if len(installed) != 0:
            while list(installed) != list(all_sections) and i < len(equipment):

                rest_equip = sorted(equipment_layers.keys(),
                                    key=lambda k: len(set(equipment_layers[k]).difference(set(installed))),
                                    reverse= True)
                rest_equip = sorted(rest_equip,
                                    key= lambda k: k[1],
                                    reverse= True)


                for re in rest_equip:
                    most_new = re
                    newlayers = [l for l in equipment_layers[most_new] if l not in installed]
                    if len(newlayers) == 0:
                        continue
                    else:
                        break

                installed.extend(newlayers)
                equips_dict[most_new[0]] = newlayers
                equips_lst.append(most_new[0])
                installed = sorted(installed)
                i += 1

            if list(installed) == list(all_sections):
                equips_lst = sorted(equips_lst)
                if equips_lst not in combinations_dict:
                    combinations_dict.append(equips_dict)
                    combinations_lst.append(equips_lst)

            all_equips_dict[f'{equips_lst}'] = equips_dict



    cost_combinations = {}
    for c in df.columns:
        if c != 'area':
            df.rename(columns={c: c.name}, inplace= True)


    for combi in combinations_dict:
        s = ''
        cost = 0
        CO2 = 0
        time = 0
        for key, value in combi.items():
            s += (key + ' + ')
            for section in value:
                d = df[key].iloc[section]
                if 'cost' in d:
                    cost += d['cost']
                if 'CO2' in d:
                    CO2 += d['CO2']
                if 'time' in d:
                    time += d['time']

        cost_combinations[s[:-3]] = {'cost': cost, 'CO2': CO2, 'time': time}

    df_optimal = pd.DataFrame.from_dict(cost_combinations, orient='index')
    df_optimal = df_optimal.sort_values(by= optimize_on)
    df_optimal['cost'] = df_optimal['cost'].round(2)
    df_optimal['CO2'] = df_optimal['CO2'].round(5)
    df_optimal['time'] = df_optimal['time'].round(5)

    return df_optimal

def combination_algorithm(dataframe, optimize_on ="cost", algorithm = 'smart_combinations', limit= None):

    if algorithm == 'smart_combinations':
        return smart_combinations(dataframe= dataframe, optimize_on= optimize_on)
    if algorithm == 'cheap_combinations' and limit == None:
        return cheap_combinations(dataframe= dataframe, optimize_on= optimize_on)
    if algorithm == 'cheap_combinations' and limit != None:
        return cheap_combinations_limit(dataframe= dataframe, optimize_on= optimize_on, limit = limit)
