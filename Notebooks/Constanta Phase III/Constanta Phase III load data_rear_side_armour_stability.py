# %%
import breakwater as bw
import pandas as pd
import os
from pathlib import Path
import numpy as np

# %%
from development_overtopping_DKA import (eurotop2018_6_5, 
                                         surf_similarity, 
                                         gamma_beta_eurotop_2018_6_9,
                                         eurotop_6_5_q, 
                                         calc_beta)
from development_armour_stability_DKA import (vangent_armour_reduction, 
                                              hudson_fixed_slope, 
                                              rock_manual_5_196, 
                                              rock_manual_5_195, 
                                              rock_manual_5_194,
                                              rock_manual_5_194_Sd,
                                              vandermeer_shallow_Sd)
from development_material_DKA import(get_class,
                                     get_Dn50)

from development_scour_DKA import(sumer_fredsoe,
                                 v_scour)
                                 
from development_toe_DKA import (toe_stability_Nod)

from breakwater.utils.exceptions import user_warning

# %%
project_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx",
                             index_col = 1,
                             sheet_name='Input_Project specific')
requirements_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx",
                                  index_col = 0,
                                  sheet_name='Input_requirements')
wave_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx",
                          index_col = 0,
                          sheet_name='input_hydrotechnical',
                         skiprows = 1)
cross_section_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx", 
                                   sheet_name='Input_Cross section',
                                  skiprows = 1)
concrete_element_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx", 
                                      sheet_name='Input_concrete_elements',
                                      index_col = 0,
                                      skiprows = 1)
gradings_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx", 
                                      sheet_name='input_rock_gradings',
                                      index_col = 0,
                                      skiprows = 2)
cross_section_validation_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx", 
                                              sheet_name='input_cross_section_validation',
                                              skiprows = 1)

# %%
#VALIDATE REAR SIDE ARMOUR STABILITY
# This calculation needs a crest height as input and therefore needs to be done after calculation of 
# the required crest height. If the calculation is done before choosing the final crest height it can give conservative
# results since the crest height in the calculation is lower than the final crest height.
# The calculation can also be performed in the verification stage only.

#To-do
# Hs-Hm0 interpretation

from numpy import cross

Calculation_case = 66
Cross_section_id = 1
# Armour = 2 layers rock, permeable core

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
tana            = cross_section_data.at[Cross_section_id, 'tan_a']
dir_structure   = cross_section_data.at[Cross_section_id, 'dir_structure']
cota_rear       = cross_section_data.at[Cross_section_id, 'cot_a_rear']
crest_roughness = cross_section_data.at[Cross_section_id, 'crest_roughness']
B               = cross_section_data.at[Cross_section_id, 'B']
z_bed           = cross_section_data.at[Cross_section_id, 'z_bed']
slope_foreshore = cross_section_data.at[Cross_section_id, 'slope_foreshore']

# Get validation data
z_crest         = cross_section_validation_data.at[Cross_section_id, 'z_crest']
grading_rear    = cross_section_validation_data.at[Cross_section_id, 'rear_armour_class']


# Get allowed Sd for the location and calculation case
LS         = wave_data.at[Calculation_case, 'Limit State']
Sd_allowed = requirements_data.at['Armour damage limit', LS]

#Intermediate calculations
h              = wl-z_bed
Rc             = z_crest - wl
beta = calc_beta(dir_structure, dir_wave)
Dn50_rear      = get_Dn50(grading_rear, gradings_data, 'Min')

waveinfo       = bw.BattjesGroenendijk(Hm0, h, slope_foreshore)
H2_per         = waveinfo.get_Hp(0.02)
Hs             = waveinfo.get_Hn(3) #NU BATTJES-GROENENDIJK VOOR Hs UIT Hmo. IS DAT WAT WE WILLEN?

#OVERRIDE wave characteristics BECAUSE W+B sheet tales Hs = Hm0
Hs = Hm0 #TO BE ADAPTED IF WE CHOOSE TO DIFFERENTIATE BETWEEN Hs AND Hm0

Delta          = (rho_a-rho_w)/rho_w
xi_s_min_1     = surf_similarity(tana, Hs, Tm_min_1, g)

# Actual calculations
Ru_1_percent, gamma_f, gamma_beta, gamma  = rock_manual_5_196(Hs, xi_s_min_1, slope_roughness, beta)
u_1_percent = rock_manual_5_195(g, slope_roughness, crest_roughness, Ru_1_percent, Rc, B, Hs)


Sd_rear = rock_manual_5_194_Sd(Dn50_rear, N, u_1_percent, Tm_min_1, Delta, cota_rear, Rc, Hs)


# print("Cross-section  : "+str(cross_section_data.at[Cross_section_id,"Structure"])+" - "+str(cross_section_data.at[Cross_section_id,"Chainage"])) # Design
# print("Wave info for  : " +str(wave_data.at[Calculation_case,"Structure"])+" - "+str(wave_data.at[Calculation_case,"Chainage"]))       # Design
# print("Limit State    : "+str(LS)) # Design
# print("Hs             : "+str(Hs)) # Design
print("xi_s_min_1     : "+str(xi_s_min_1))
print("slope_roughness: "+str(slope_roughness))
print("crest_roughness: "+str(crest_roughness))
print("beta           : "+str(beta))
print("gamma_f        : "+str(gamma_f))
print("gamma_beta     : "+str(gamma_beta))
print("gamma          : "+str(gamma))
print("Rc             : "+str(Rc))
print("B              : "+str(B))
print("Sd_allowed     : "+str(Sd_allowed))
print("N              : "+str(N))
print("Tm_min_1       : "+str(Tm_min_1))
print("Delta          : "+str(Delta))
print("cota_rear:     : "+str(cota_rear))
print("u_1_percent    : "+str(u_1_percent)) # Design
print("Ru_1_percent   : "+str(Ru_1_percent)) # Design
print("Dn50_rear      : "+str(Dn50_rear)) # Design

