# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:48:09 2022

@author: JY6
"""
import pandas as pd
import numpy as np
from breakwater.utils.exceptions import user_warning
# gradings_data = pd.read_excel(r"C:\Users\JY6\github\breakwater\Notebooks\Constanta Phase III\Input data\test_data_phase_II.xlsx", 
#                                       sheet_name='input_rock_gradings',
#                                       index_col = 0,
#                                       skiprows = 2)


# Dn50 = 0.6
# rho_a = 2600


# M_required = rho_a * Dn50 ** 3


def get_class(Dn50, rho_a, gradings_table):
    """Get the rock class for a given Dn50

    Parameters
    ----------
    Dn50 : float
        nominal diameter of the armourstone [m]
    rho_a : float
        armourstone density [kg/m3]
    gradings_table : DataFrame
        table with the armourstone properties, based on Constanta Phase III
        input file

    Returns
    -------
    str
        Rock class

    Raises
    ------
    RockGradingError
        If the computed Dn50 of the armour layer is out of range for
        the specified rock grading
    """
    rock_class = None
    M_required = rho_a * Dn50 ** 3
    
    # rock_class = gradings_data[gradings_data['W_M50']>M_required].iloc[0].name
    
    for index, row in gradings_table.iterrows():
        if row['W_M50']>= M_required:
            rock_class = index
            break
        
    if rock_class == None:
        max_class = gradings_table['W_M50'].index[-1]
        #max_mass  = gradings_table[max_class]["W_M50"][1]
        max_dn = np.round(gradings_table.at[max_class,'Dn50 min'], 3)
        rock_class = max_class
        user_warning(
                f"Dn50 = {np.round(Dn50, 3)} is out of range for the specified"
                f" rock grading, {max_dn} m is the maximum possible Dn50. {rock_class} is taken as "
                f" the grading but does not apply to the filter rules."
            )
    return rock_class

def get_Dn50(rock_class, gradings_table):
    """Get the rock class for a given Dn50

    Parameters
    ----------
    rock_class : str
        class that the average Dn50 is selected for
    gradings_table : DataFrame
        table with the armourstone properties, based on Constanta Phase III
        input file

    Returns
    -------
    str
        Rock class

    Raises
    ------
    RockGradingError
        If the computed Dn50 of the armour layer is out of range for
        the specified rock grading
    """
    
    Dn50 = gradings_table.at[rock_class, 'Dn50 av']
    
    return Dn50

# rock_class = get_class(Dn50, rho_a, gradings_data)
# Dn50_class = get_Dn50(rock_class, gradings_data)
