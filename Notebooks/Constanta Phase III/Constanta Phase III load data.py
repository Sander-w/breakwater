# %%
import breakwater as bw
import pandas as pd
import os
from pathlib import Path
import numpy as np

# %%
from development_overtopping_DKA import eurotop2018_6_5, surf_similarity, gamma_beta_eurotop_2018_6_9


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


# armour_layer = 'Rock'
# layers       = 2
# permeability = 'permeable'

# wave_data["xi_m_min_1"] = wave_data.apply(
#     lambda x: surf_similarity(
#         get_cross_section_data(x["Location"])[0], 
#         x["Hm0"], 
#         x["Tm-1,0"], 
#         g
#     ), 
#     axis=1
# )

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
            beta       = dir_wave-(dir_structure-360) #Done like this for ECn2 north. Should find more robust solution
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

wave_data.to_excel("wave_data_intermediate.xlsx")

results = []
for chainage in wave_data.Chainage.unique():
    chainage_data = []
    
    chainage_data.append(chainage)

    max_LS = max(wave_data[wave_data["Chainage"] == chainage]["Rock_public access_z_crest"])
    chainage_data.append(list(wave_data[wave_data["Rock_public access_z_crest"] == max_LS]["Structure"])[0])
    chainage_data.append(list(wave_data[wave_data["Rock_public access_z_crest"] == max_LS]["Limit State"])[0])
    chainage_data.append(list(wave_data[wave_data["Rock_public access_z_crest"] == max_LS]["Offshore bin"])[0])
    chainage_data.append(list(wave_data[wave_data["Rock_public access_z_crest"] == max_LS]["Hm0"])[0])
    chainage_data.append(max_LS)

    max_LS = max(wave_data[wave_data["Chainage"] == chainage]["Rock_restricted access_z_crest"])
    chainage_data.append(list(wave_data[wave_data["Rock_restricted access_z_crest"] == max_LS]["Limit State"])[0])
    chainage_data.append(list(wave_data[wave_data["Rock_restricted access_z_crest"] == max_LS]["Offshore bin"])[0])
    chainage_data.append(list(wave_data[wave_data["Rock_restricted access_z_crest"] == max_LS]["Hm0"])[0])
    chainage_data.append(max_LS)

    max_LS = max(wave_data[wave_data["Chainage"] == chainage]["Xbloc_public access_z_crest"])
    chainage_data.append(list(wave_data[wave_data["Xbloc_public access_z_crest"] == max_LS]["Limit State"])[0])
    chainage_data.append(list(wave_data[wave_data["Xbloc_public access_z_crest"] == max_LS]["Offshore bin"])[0])
    chainage_data.append(list(wave_data[wave_data["Xbloc_public access_z_crest"] == max_LS]["Hm0"])[0])
    chainage_data.append(max_LS)

    max_LS = max(wave_data[wave_data["Chainage"] == chainage]["Xbloc_restricted access_z_crest"])
    chainage_data.append(list(wave_data[wave_data["Xbloc_restricted access_z_crest"] == max_LS]["Limit State"])[0])
    chainage_data.append(list(wave_data[wave_data["Xbloc_restricted access_z_crest"] == max_LS]["Offshore bin"])[0])
    chainage_data.append(list(wave_data[wave_data["Xbloc_restricted access_z_crest"] == max_LS]["Hm0"])[0])
    chainage_data.append(max_LS)

    results.append(chainage_data)

print(results)
columns = [
    "Chainage",
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
results_df.to_excel("wave_data_z_crest_1.xlsx")




# %%
# CALCULATE REQUIRED STONE DIAMETER (Rock armour)

# Calculation_case = 57
# Cross_section_id = 0
# # Armour = 2 layers rock, permeable core

# Storm_duration = 6 #Hours
# Safety = 0 #Safety factor on top of regular VDM constants

# # Open project specific parameters
# g     = project_data.at['g', 'Value']
# rho_a = project_data.at['rho_r', 'Value'] #For rock armour. Needs adaptation for concrete
# rho_w = project_data.at['rho_w', 'Value']
# Cpl_shallow = project_data.at['Cpl_shallow', 'Value']
# Cs_shallow = project_data.at['Cs_shallow', 'Value']

# #Get info for sea state
# Hm0      = wave_data.at[Calculation_case, 'Hm0']
# wl       = wave_data.at[Calculation_case, 'wl']
# Tm_min_1 = wave_data.at[Calculation_case, 'Tm-1,0']
# dir_wave = wave_data.at[Calculation_case, 'dir_wave']
# Tm_0_2   = wave_data.at[Calculation_case, 'Tm0,2']
# N        = np.round(Storm_duration*3600/Tm_0_2)


# # Open structure specific parameters
# tana            = cross_section_data.at[Cross_section_id, 'tan_a']
# dir_structure   = cross_section_data.at[Cross_section_id, 'dir_structure']
# safety          = cross_section_data.at[Cross_section_id, 'safety']
# z_bed           = cross_section_data.at[Cross_section_id, 'z_bed']
# slope_foreshore = cross_section_data.at[Cross_section_id, 'slope_foreshore']
# P               = cross_section_data.at[Cross_section_id, 'P']


# # Get allowed Sd for the location and calculation case
# LS         = wave_data.at[Calculation_case, 'Limit State']
# Sd_allowed = requirements_data.at['Armour damage limit', LS]

# #Intermediate calculations
# h              = wl-z_bed
# waveinfo       = bw.BattjesGroenendijk(Hm0, h, slope_foreshore)
# H2_per         = waveinfo.get_Hp(0.02)
# Hs             = waveinfo.get_Hn(3) #NU BATTJES-GROENENDIJK VOOR Hs UIT Hmo. IS DAT WAT WE WILLEN?

# #OVERRIDE wave characteristics BECAUSE THERE IS A DIFFERENCE BETWEEN W+B SHEET AND BW TOOL
# Hs = Hm0
# H2_per = Hs*1.34
# print("Warning: Hs and H2_per taken from W+B calculation")

# Delta          = (rho_a-rho_w)/rho_w
# xi_s_min_1     = surf_similarity(tana, Hs, Tm_0_2, g)
# alpha          = np.arctan(tana)



# # Check validity of Van der Meer shallow
# vdm_shallow_validity = h/Hs

# #Calculate required stone diameter without reduction
# Dn50 = bw.core.vandermeer_shallow(Hs, 
#                                   H2_per, 
#                                   Delta, 
#                                   P, 
#                                   Sd_allowed, 
#                                   N, 
#                                   xi_s_min_1, 
#                                   alpha, 
#                                   Cpl = Cpl_shallow, 
#                                   Cs = Cs_shallow, 
#                                   safety = Safety)

# #Reduce Dn50 with reduction factors from DAR


# # %%
# print("Cross-section : "+str(cross_section_data.at[Cross_section_id,"Structure"])+" - "+str(cross_section_data.at[Cross_section_id,"Chainage"]))
# print("Wave info for : " +str(wave_data.at[Calculation_case,"Structure"])+" - "+str(wave_data.at[Calculation_case,"Chainage"]))      
# print("Limit State   : "+str(LS))
# print("Allowed Sd    : "+str(Sd_allowed))
# print("g         : "+str(g))
# print("h         : "+str(h))
# print("Delta     : "+str(Delta))
# print("P         : "+str(P))
# print("Tangent a : "+str(tana))
# print("Hs        : "+str(Hs))
# print("H2%       : "+str(H2_per))
# print("H2%/Hs    : "+str(H2_per/Hs))
# print("Tm0,2     : "+str(Tm_0_2))
# print("N         : "+str(N))
# print("xi_s_min_1: "+str(xi_s_min_1))
# print("Dn50      : "+str(Dn50))
# print("h/Hs (<3) : "+str(vdm_shallow_validity))


# # %%
# Hs/Hm0


# # %%


# # %%


# # %%


# # %%



