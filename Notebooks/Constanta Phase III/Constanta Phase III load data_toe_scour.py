#%% Define input and output files
input_file = "test_data_phase_II.xlsx"
output_file = 'An1_first_test_results.xlsx'

# %% Import functions and packages
import breakwater as bw
import pandas as pd
import os
from pathlib import Path
import numpy as np
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)
logging.info("Initiated script")

# %% Import local functions
from development_overtopping_DKA import (eurotop2018_6_5,
                                         surf_similarity,
                                         gamma_beta_eurotop_2018_6_9,
                                         calc_beta)
from development_armour_stability_DKA import (vangent_armour_reduction,
                                              hudson_fixed_slope,
                                              rock_manual_5_196,
                                              rock_manual_5_195,
                                              rock_manual_5_194)
from development_material_DKA import(get_class,
                                     get_Dn50)

from development_scour_DKA import(sumer_fredsoe,
                                 v_scour)


from breakwater.utils.exceptions import user_warning

# %% Import input data
project_data = pd.read_excel(Path("./Input data/") / input_file,
                             index_col = 1,
                             sheet_name='Input_Project specific')
requirements_data = pd.read_excel(Path("./Input data/") / input_file,
                                  index_col = 0,
                                  sheet_name='Input_requirements')

wave_data = pd.read_excel(Path("./Input data/") / input_file,
                          index_col = 0,
                          sheet_name='input_hydrotechnical',
                         skiprows = 1)
wave_data["Location"] = wave_data["Structure"] + wave_data["Chainage"]
columns = wave_data.columns.tolist()[:-1]
columns.insert(2,"Location")
wave_data = wave_data[columns]

cross_section_data = pd.read_excel(Path("./Input data/") / input_file,
                                   sheet_name='Input_Cross section',
                                  skiprows = 1)
cross_section_data["Location"] = cross_section_data["Structure"] + cross_section_data["Chainage"]
cross_section_data = cross_section_data.set_index('Location')
concrete_element_data = pd.read_excel(Path("./Input data/") / input_file,
                                      sheet_name='Input_concrete_elements',
                                      index_col = 0,
                                      skiprows = 1)
gradings_data = pd.read_excel(Path("./Input data/") / input_file,
                                      sheet_name='input_rock_gradings',
                                      index_col = 0,
                                      skiprows = 2)

# %% CALCULATE TOE STABILITY
[].insert

Nod_allowed_list = []
Delta_list = []
Hs_list = []
h_list = []
t_filter_list = []
toe_layer_thickness_list = []
ht_list = []
Dn50_toe_list = []
toe_armour_class_list = []
Dn50_toe_temp_list = []

for Calculation_case in range(1, len(wave_data.index)+1):
    # Calculation_case = 1
    Cross_section_id = wave_data.at[Calculation_case, 'Location']

    # Open project specific parameters
    g            = project_data.at['g'           , 'Value']
    rho_a        = project_data.at['rho_r'       , 'Value']
    rho_w        = project_data.at['rho_w'       , 'Value']
    Lt           = project_data.at['Lt'          , 'Value']
    t_underlayer = project_data.at['t_underlayer', 'Value']

    #Get info for sea state
    Hm0      = wave_data.at[Calculation_case, 'Hm0']
    wl       = wave_data.at[Calculation_case, 'wl']

    # Open structure specific parameters
    z_bed           = cross_section_data.at[Cross_section_id, 'z_bed']

    # Get allowed Sd for the location and calculation case
    LS          = wave_data.at[Calculation_case, 'Limit State']
    Nod_allowed = requirements_data.at['Toe damage limit', LS]

    # Intermediate parameters
    h               = wl-z_bed
    Delta          = (rho_a-rho_w)/rho_w

    #TEMPORARY FIX TO BE IN LINE WITH PHASE II
    Hs = Hm0

    # Calculate required stone size for toe

    # use while loop since ht depends on Dn50 of the toe
    # first set temporary values for the while loop
    Dn50_toe = 0
    Dn50_toe_temp = 0
    compute_toe = True
    counter = 0
    t_filter = 0

    toe_layer_thickness = t_underlayer + t_filter + Dn50_toe*Lt
    ht = h - toe_layer_thickness

    while compute_toe:
        Dn50_toe_computed = bw.core.toe_stability(Hs,
                                                h,
                                                ht,
                                                Delta,
                                                Nod_allowed)
        if np.isnan(Dn50_toe_computed):
            Dn50_toe_computed = 0.01

        # check for convergence
        if Dn50_toe_computed - Dn50_toe_temp < 0 or counter > 50:
            # value has converged, so break loop
            compute_toe = False

        counter += 1

        # check if computed Dn50 is larger than current normative Dn50
        if Dn50_toe_computed > Dn50_toe:
            # if larger the normative Dn50 must be changed
            Dn50_toe = Dn50_toe_computed

        # Set Dn50 toe class average
        class_toe = get_class(Dn50_toe_computed, rho_a, gradings_data)
        Dn50_toe_computed = get_Dn50(class_toe, gradings_data, 'Av')
        # Check if filter layer needs to be applied
        if gradings_data.at[class_toe, 'Toe underlayer'] in gradings_data.index:
            class_filter = gradings_data.at[class_toe, 'Toe underlayer']
            Dn50_filter = get_Dn50(class_filter, gradings_data, 'Av')
            t_filter = Dn50_filter*Lt




        # replace old values with the new ones
        Dn50_toe_temp = Dn50_toe_computed
        toe_layer_thickness = t_underlayer + t_filter + Dn50_toe_temp*Lt
        ht = h - toe_layer_thickness

    toe_armour_class = get_class(Dn50_toe, rho_a, gradings_data)

    Nod_allowed_list.append(Nod_allowed)
    Delta_list.append(Delta)
    Hs_list.append(Hs)
    h_list.append(h)
    t_filter_list.append(t_filter)
    toe_layer_thickness_list.append(toe_layer_thickness)
    ht_list.append(ht)
    Dn50_toe_list.append(Dn50_toe)
    toe_armour_class_list.append(toe_armour_class)
    Dn50_toe_temp_list.append(Dn50_toe_temp)

wave_data["Nod_allowed"]= Nod_allowed_list
wave_data["Delta"]= Delta_list
wave_data["Hs"]= Hs_list
wave_data["h"]= h_list
wave_data["t_filter"]= t_filter_list
wave_data["toe_layer_thickness"]= toe_layer_thickness_list
wave_data["ht"]= ht_list
wave_data["Dn50_toe_calculated"]= Dn50_toe_list
wave_data["toe_armour_class"]= toe_armour_class_list
wave_data["Dn50_grading"]= Dn50_toe_temp_list

wave_data.to_excel("wave_data_intermediate_toe_stability.xlsx")

logging.info("Finished toe stability intermediate section")


results = []
for location in wave_data.Location.unique():
    location_summary = []

    location_summary.append(location)

    normative_case = max(wave_data[wave_data["Location"] == location]["Dn50_grading"])
    location_summary.append(list(wave_data[wave_data["Dn50_grading"] == normative_case]["Structure"])[0])
    location_summary.append(list(wave_data[wave_data["Dn50_grading"] == normative_case]["Limit State"])[0])
    location_summary.append(list(wave_data[wave_data["Dn50_grading"] == normative_case]["Offshore bin"])[0])
    location_summary.append(normative_case)

    results.append(location_summary)

columns = [
    "Location",
    "Structure",
    "LS",
    "Offshore bin",
    "max Dn50_grading",
]
results_df = pd.DataFrame(results, columns=columns)
# results_df.to_excel("wave_data_design_toe_stability.xlsx")

logging.info("Finished toe stability design section")

# %% CALCULATE SCOUR DEPTH

wave_data = pd.read_excel(Path("./Input data/") / input_file,
                          index_col = 0,
                          sheet_name='input_hydrotechnical',
                         skiprows = 1)
wave_data["Location"] = wave_data["Structure"] + wave_data["Chainage"]
columns = wave_data.columns.tolist()[:-1]
columns.insert(2,"Location")
wave_data = wave_data[columns]

Hs_list = []
Tp_list = []
h_list = []
S_list = []

for Calculation_case in range(1, len(wave_data.index)+1):
    # Calculation_case = 1
    Cross_section_id = wave_data.at[Calculation_case, 'Location']


    # Open project specific parameters
    g            = project_data.at['g'           , 'Value']
    C2_sf        = project_data.at['C2_sf'       , 'Value']

    #Get info for sea state
    Hm0      = wave_data.at[Calculation_case, 'Hm0']
    Tp       = wave_data.at[Calculation_case, 'Tp']
    wl       = wave_data.at[Calculation_case, 'wl']

    # Open structure specific parameters
    z_bed           = cross_section_data.at[Cross_section_id, 'z_bed']

    # Intermediate parameters
    h               = wl-z_bed

    #TEMPORARY FIX TO BE IN LINE WITH PHASE II
    Hs = Hm0

    # Calculate scour depth
    S = sumer_fredsoe(Hs, Tp, g, h, C2_sf)

    Hs_list.append(Hs)
    Tp_list.append(Tp)
    h_list.append(h)
    S_list.append(S)

wave_data['Hs'] = Hs_list
wave_data['Tp'] = Tp_list
wave_data['h'] = h_list
wave_data['S'] = S_list

wave_data.to_excel("wave_data_intermediate_scour.xlsx")

logging.info("Finished scour intermediate section")


results = []
for location in wave_data.Location.unique():
    location_summary = []

    location_summary.append(location)

    normative_case = max(wave_data[wave_data["Location"] == location]["S"])
    location_summary.append(list(wave_data[wave_data["S"] == normative_case]["Structure"])[0])
    location_summary.append(list(wave_data[wave_data["S"] == normative_case]["Limit State"])[0])
    location_summary.append(list(wave_data[wave_data["S"] == normative_case]["Offshore bin"])[0])
    location_summary.append(normative_case)

    results.append(location_summary)

columns = [
    "Location",
    "Structure",
    "LS",
    "Offshore bin",
    "max S",
]
results_df = pd.DataFrame(results, columns=columns)
results_df.to_excel("wave_data_design_scour.xlsx")

logging.info("Finished scour design section")

# #%% Write to single excel file per structure
# # If the file exists: load file and adapt tabs only. If the file does not exist: create it

# if os.path.exists(output_file):
#     writer = pd.ExcelWriter(output_file,
#                             engine = 'openpyxl',
#                             mode = 'a',
#                             if_sheet_exists = 'replace')
# else:
#     writer = pd.ExcelWriter(output_file,
#                             engine = 'openpyxl',
#                             mode = 'w')

# wave_data.to_excel(writer, sheet_name = 'armour_sta_intermediate', index = False)
# results_df.to_excel(writer, sheet_name = 'armour_sta_summary', index = False)

# #writer.save()
# writer.close()




# %%
