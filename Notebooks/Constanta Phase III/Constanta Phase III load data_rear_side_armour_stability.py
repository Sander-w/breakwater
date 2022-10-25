#%% Define input and output files
project_parameters = "project_specific_input.xlsx"
input_file = "first_test_data_phase_III.xlsx"
output_file = 'An1_first_test_results.xlsx'

# %% Import functions and packages
import breakwater as bw
import pandas as pd
import os
from pathlib import Path
import numpy as np
from openpyxl import load_workbook
from os import path

# %%
from development_overtopping_DKA import eurotop2018_6_5, surf_similarity, gamma_beta_eurotop_2018_6_9, calc_beta
from development_armour_stability_DKA import rock_manual_5_194, rock_manual_5_194_Sd, rock_manual_5_195, rock_manual_5_196
from development_material_DKA import get_class

# %% Import input data
project_data = pd.read_excel(Path("./Input data/") / project_parameters,
                             index_col = 1,
                             sheet_name='Input_Project specific')
requirements_data = pd.read_excel(Path("./Input data/") / project_parameters,
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
gradings_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx", 
                                      sheet_name='input_rock_gradings',
                                      index_col = 0,
                                      skiprows = 2)


# %%
#VALIDATE REAR SIDE ARMOUR STABILITY
# This calculation needs a crest height as input and therefore needs to be done after calculation of 
# the required crest height. If the calculation is done before choosing the final crest height it can give conservative
# results since the crest height in the calculation is lower than the final crest height.
# The calculation can also be performed in the verification stage only.

#To-do
# Hs-Hm0 interpretation

xi_s_min_1_list = []
slope_roughness_list = []
crest_roughness_list = []
beta_list = []
gamma_f_list = []
gamma_beta_list = []
gamma_list = []
Rc_list = []
B_list = []
Sd_allowed_list = []
N_list = []
Tm_min_1_list = []
Delta_list = []
cota_rear_list = []
u_1_percent_list = []
Ru_1_percent_list = []
Dn50_rear_list = []


for Calculation_case in range(1, len(wave_data.index)+1):
    # Calculation_case = 7
    
    Cross_section_id = wave_data.at[Calculation_case, 'Location']

    # Open project specific parameters
    g               = project_data.at['g'              , 'Value']
    rho_a           = project_data.at['rho_r'          , 'Value']
    rho_w           = project_data.at['rho_w'          , 'Value']
    slope_roughness = project_data.at['slope_roughness', 'Value']
    Storm_duration  = project_data.at['storm_duration' , 'Value']

    #Get info for sea state
    Hm0      = wave_data.at[Calculation_case, 'Hm0']
    dir_wave = wave_data.at[Calculation_case, 'dir_wave']
    Tm_0_2   = wave_data.at[Calculation_case, 'Tm0,2']
    Tm_min_1 = wave_data.at[Calculation_case, 'Tm-1,0']
    wl       = wave_data.at[Calculation_case, 'wl']
    N        = np.round(Storm_duration*3600/Tm_0_2)

    # Open structure specific parameters
    tana            = cross_section_data.at[Cross_section_id, 'tan_a_rock']
    dir_structure   = cross_section_data.at[Cross_section_id, 'dir_structure']
    cota_rear       = cross_section_data.at[Cross_section_id, 'cot_a_rear']
    crest_roughness = cross_section_data.at[Cross_section_id, 'crest_roughness']
    B               = cross_section_data.at[Cross_section_id, 'B']
    z_bed           = cross_section_data.at[Cross_section_id, 'z_bed']
    slope_foreshore = cross_section_data.at[Cross_section_id, 'slope_foreshore']

    # Get allowed Sd for the location and calculation case
    LS         = wave_data.at[Calculation_case, 'Limit State']
    Sd_allowed = requirements_data.at['Armour damage limit', LS]

    #Intermediate calculations
    h              = wl-z_bed
    beta = calc_beta(dir_structure, dir_wave)
    waveinfo       = bw.BattjesGroenendijk(Hm0, h, slope_foreshore)
    H2_per         = waveinfo.get_Hp(0.02)
    Hs             = waveinfo.get_Hn(3) #NU BATTJES-GROENENDIJK VOOR Hs UIT Hmo. IS DAT WAT WE WILLEN?

    #OVERRIDE wave characteristics BECAUSE W+B sheet tales Hs = Hm0
    Hs = Hm0 #TO BE ADAPTED IF WE CHOOSE TO DIFFERENTIATE BETWEEN Hs AND Hm0

    Delta          = (rho_a-rho_w)/rho_w
    xi_s_min_1     = surf_similarity(tana, Hs, Tm_min_1, g)

    Rc = 1.99-wl

    # Actual calculations
    Ru_1_percent, gamma_f, gamma_beta, gamma  = rock_manual_5_196(Hs, xi_s_min_1, slope_roughness, beta)
    u_1_percent = rock_manual_5_195(g, slope_roughness, crest_roughness, Ru_1_percent, Rc, B, Hs)
    #Rc IS NOW JUST THE PRESENT Rc. THIS NEEDS TO BE COUPLED TO CALCULATION CASE
    Dn50_rear = rock_manual_5_194(Sd_allowed, N, u_1_percent, Tm_min_1, Delta, cota_rear, Rc, Hs)

    # Select grading    
    grading_rear = get_class(Dn50_rear, rho_a, gradings_data)

    xi_s_min_1_list.append(xi_s_min_1)
    slope_roughness_list.append(slope_roughness)
    crest_roughness_list.append(crest_roughness)
    beta_list.append(beta)
    gamma_f_list.append(gamma_f)
    gamma_beta_list.append(gamma_beta)
    gamma_list.append(gamma)
    Rc_list.append(Rc)
    B_list.append(B)
    Sd_allowed_list.append(Sd_allowed)
    N_list.append(N)
    Tm_min_1_list.append(Tm_min_1)
    Delta_list.append(Delta)
    cota_rear_list.append(cota_rear)
    u_1_percent_list.append(u_1_percent)
    Ru_1_percent_list.append(Ru_1_percent)
    Dn50_rear_list.append(Dn50_rear)

wave_data['xi_s_min_1'] = xi_s_min_1_list
wave_data['slope_roughness'] = slope_roughness_list
wave_data['crest_roughness'] = crest_roughness_list
wave_data['beta'] = beta_list
wave_data['gamma_f'] = gamma_f_list
wave_data['gamma_beta'] = gamma_beta_list
wave_data['gamma'] = gamma_list
wave_data['Rc'] = Rc_list
wave_data['B'] = B_list
wave_data['Sd_allowed'] = Sd_allowed_list
wave_data['N'] = N_list
wave_data['Tm_min_1'] = Tm_min_1_list
wave_data['Delta'] = Delta_list
wave_data['cota_rear'] = cota_rear_list
wave_data['u_1_percent'] = u_1_percent_list
wave_data['Ru_1_percent'] = Ru_1_percent_list
wave_data['Dn50_rear'] = Dn50_rear_list

wave_data.to_excel("wave_data_intermediate_rear_side_armour_stability.xlsx")

results = []
for location in wave_data.Location.unique():
    location_summary = []
    
    location_summary.append(location)

    print(location)

    normative_case = max(wave_data[wave_data["Location"] == location]["Dn50_rear"])
    location_summary.append(list(wave_data[wave_data["Dn50_rear"] == normative_case]["Structure"])[0])
    location_summary.append(list(wave_data[wave_data["Dn50_rear"] == normative_case]["Limit State"])[0])
    location_summary.append(list(wave_data[wave_data["Dn50_rear"] == normative_case]["Offshore bin"])[0])
    location_summary.append(list(wave_data[wave_data["Dn50_rear"] == normative_case]["u_1_percent"])[0])
    location_summary.append(list(wave_data[wave_data["Dn50_rear"] == normative_case]["Ru_1_percent"])[0])
    location_summary.append(normative_case)

    results.append(location_summary)

columns = [
    "Location",
    "Structure",
    "LS",
    "Offshore bin",
    "u_1_percent",
    "Ru_1_percent",
    "max Dn50_rear",
]
results_df = pd.DataFrame(results, columns=columns)
# results_df.to_excel("wave_data_design_toe_stability.xlsx")


if os.path.exists(output_file):
    writer = pd.ExcelWriter(output_file,
                            engine = 'openpyxl',
                            mode = 'a',
                            if_sheet_exists = 'replace')
else:
    writer = pd.ExcelWriter(output_file,
                            engine = 'openpyxl',
                            mode = 'w')

wave_data.to_excel(writer, sheet_name = 'rear_side_armour_stability_intermediate', index = False)
results_df.to_excel(writer, sheet_name = 'rear_side_armour_stability_summary', index = False)

#writer.save()
writer.close()
# %%
