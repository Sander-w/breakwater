# %%
import breakwater as bw
import pandas as pd
import os
from pathlib import Path
import numpy as np

# %%
from development_overtopping_DKA import eurotop2018_6_5, surf_similarity, gamma_beta_eurotop_2018_6_9, calc_beta


# %%
project_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx",
    index_col = 1,
    sheet_name='Input_Project specific')
requirements_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx",
    index_col = 0,
    sheet_name='Input_requirements')
wave_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx",
    # index_col = 0,
    sheet_name='input_hydrotechnical',
    skiprows = 1)
wave_data["Location"] = wave_data["Structure"] + wave_data["Chainage"]
columns = wave_data.columns.tolist()[:-1]
columns.insert(2,"Location")
wave_data = wave_data[columns]
cross_section_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx", 
    sheet_name='Input_Cross section',
    skiprows = 1)
cross_section_data["Location"] = cross_section_data["Structure"] + cross_section_data["Chainage"]
cross_section_data = cross_section_data.set_index('Location')


# %%
# CALCULATE CREST HEIGHT FOR OVERTOPPING

# Do we want to put this entire part in a function to keep the scripts clean a bit?

# Still necessary: Loop over all calculation cases for a cross section

# Still necessary: Error/notification if there is no overtopping limit

# Still necessary: implementation to select correct overtopping limit depending
# on if public access is allowed

# Still necessary: implementation to select correct roughness
# depending on armour case.

# Calculation_case = 3
   

def get_cross_section_data(location):
    tana          = cross_section_data.at[location, 'tan_a']
    dir_structure = cross_section_data.at[location, 'dir_structure']
    safety        = cross_section_data.at[location, 'safety']
    return tana, dir_structure, safety

def get_requirements_data(access, LS):
    return requirements_data.at[access, LS]

g = project_data.at['g', 'Value']

for armour_layer in ["Rock", "Xbloc"]:
    for access in ["Overtopping limit public access", "Overtopping limit restricted access"]:
        q_allowed_list = []
        xi_m_min_1_list = []
        Rc_list = []
        beta_list = []
        gamma_f_list = []
        gamma_beta_list = []
        z_crest_list = []
    
        for Calculation_case in range(0, len(wave_data.Calculation_case)):
            Cross_section_id = wave_data.at[Calculation_case, 'Location']
            # Public access: allowed
            # Armour = 2 layers rock, permeable core

            # Open project specific parameters
            

            # Open structure specific parameters
            tana, dir_structure, safety = get_cross_section_data(Cross_section_id)

            # Get info for sea state. Currently implemented to do only the first sea state
            Hm0      = wave_data.at[Calculation_case, 'Hm0']
            wl       = wave_data.at[Calculation_case, 'wl']
            Tm_min_1 = wave_data.at[Calculation_case, 'Tm-1,0']
            dir_wave = wave_data.at[Calculation_case, 'dir_wave']

            # Get allowed q for the location and calculation case
            LS        = wave_data.at[Calculation_case, 'Limit State']
            q_allowed = get_requirements_data(access, LS)


            #Calculate roughness reduction
            # Set armour for rougness reduction

            # Calculate surf-similarity with Hm0 and Tm-1,0
            xi_m_min_1 = surf_similarity(tana, Hm0, Tm_min_1, g)

            
            if armour_layer == "Xbloc":
                gamma_f = bw.core.overtopping.gamma_f(armour_layer, xi_m_min_1)
            elif armour_layer == "Rock":
                gamma_f = bw.core.overtopping.gamma_f(armour_layer, xi_m_min_1, layers = '2', permeability = 'permeable')
            #Calculate obliqueness reduction
            beta       = calc_beta(dir_structure, dir_wave)
            gamma_beta = gamma_beta_eurotop_2018_6_9(beta)


            Rc = eurotop2018_6_5(g, Hm0, q_allowed, gamma_f, gamma_beta, limit = False, safety=safety)
            z_crest = Rc + wl

            q_allowed_list.append(q_allowed)
            xi_m_min_1_list.append(xi_m_min_1)
            Rc_list.append(Rc)
            beta_list.append(beta)
            gamma_f_list.append(gamma_f)
            gamma_beta_list.append(gamma_beta)
            z_crest_list.append(z_crest)

        combined_string = armour_layer + "_" + access.split("Overtopping limit ")[1]        
        wave_data["xi_m_min_1"] = xi_m_min_1_list        
        wave_data["beta"] = beta_list        
        wave_data["gamma_beta"] = gamma_beta_list
        wave_data[combined_string + "_q_allowed"] = q_allowed_list
        wave_data[combined_string + "_gamma_f"] = gamma_f_list
        wave_data[combined_string + "_Rc"] = Rc_list
        wave_data[combined_string + "_z_crest"] = z_crest_list

wave_data.to_excel("wave_data_intermediate_crest_height.xlsx")

results = []
for location in wave_data.Location.unique():
    location_summary = []
    
    location_summary.append(location)

    normative_case = max(wave_data[wave_data["Location"] == location]["Rock_public access_z_crest"])
    location_summary.append(list(wave_data[wave_data["Rock_public access_z_crest"] == normative_case]["Structure"])[0])
    location_summary.append(list(wave_data[wave_data["Rock_public access_z_crest"] == normative_case]["Limit State"])[0])
    location_summary.append(list(wave_data[wave_data["Rock_public access_z_crest"] == normative_case]["Offshore bin"])[0])
    location_summary.append(list(wave_data[wave_data["Rock_public access_z_crest"] == normative_case]["Hm0"])[0])
    location_summary.append(normative_case)

    normative_case = max(wave_data[wave_data["Location"] == location]["Rock_restricted access_z_crest"])
    location_summary.append(list(wave_data[wave_data["Rock_restricted access_z_crest"] == normative_case]["Limit State"])[0])
    location_summary.append(list(wave_data[wave_data["Rock_restricted access_z_crest"] == normative_case]["Offshore bin"])[0])
    location_summary.append(list(wave_data[wave_data["Rock_restricted access_z_crest"] == normative_case]["Hm0"])[0])
    location_summary.append(normative_case)

    normative_case = max(wave_data[wave_data["Location"] == location]["Xbloc_public access_z_crest"])
    location_summary.append(list(wave_data[wave_data["Xbloc_public access_z_crest"] == normative_case]["Limit State"])[0])
    location_summary.append(list(wave_data[wave_data["Xbloc_public access_z_crest"] == normative_case]["Offshore bin"])[0])
    location_summary.append(list(wave_data[wave_data["Xbloc_public access_z_crest"] == normative_case]["Hm0"])[0])
    location_summary.append(normative_case)

    normative_case = max(wave_data[wave_data["Location"] == location]["Xbloc_restricted access_z_crest"])
    location_summary.append(list(wave_data[wave_data["Xbloc_restricted access_z_crest"] == normative_case]["Limit State"])[0])
    location_summary.append(list(wave_data[wave_data["Xbloc_restricted access_z_crest"] == normative_case]["Offshore bin"])[0])
    location_summary.append(list(wave_data[wave_data["Xbloc_restricted access_z_crest"] == normative_case]["Hm0"])[0])
    location_summary.append(normative_case)

    results.append(location_summary)

print(results)
columns = [
    "Location",
    "Structure",
    "Rock, Public, LS", 
    "Rock, Public, Offshore bin", 
    "Rock, Public, Hm0", 
    "Rock, Public, max z_crest", 
    "Rock, Restricted, LS", 
    "Rock, Restricted, Offshore bin",
    "Rock, Restricted, Hm0", 
    "Rock, Restricted, max z_crest", 
    "Xbloc, Public, LS", 
    "Xbloc, Public, Offshore", 
    "Xbloc, Public, Hm0", 
    "Xbloc, Public, max z_crest", 
    "Xbloc, Restricted, LS", 
    "Xbloc, Restricted, Offshore", 
    "Xbloc, Restricted, Hm0", 
    "Xbloc, Restricted, max z_crest"
]
results_df = pd.DataFrame(results, columns=columns)
results_df.to_excel("wave_data_design_crest_height.xlsx")