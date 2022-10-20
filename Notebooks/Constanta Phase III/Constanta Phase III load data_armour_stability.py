# %%
import breakwater as bw
import pandas as pd
import os
from pathlib import Path
import numpy as np

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)
logging.info("Initiated script")

# %%
from development_overtopping_DKA import (eurotop2018_6_5, 
    surf_similarity,
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
wave_data["Location"] = wave_data["Structure"] + wave_data["Chainage"]
columns = wave_data.columns.tolist()[:-1]
columns.insert(2,"Location")
wave_data = wave_data[columns]
cross_section_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx", 
    sheet_name='Input_Cross section',
    skiprows = 1)
cross_section_data["Location"] = cross_section_data["Structure"] + cross_section_data["Chainage"]
cross_section_data = cross_section_data.set_index('Location')
concrete_element_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx", 
    sheet_name='Input_concrete_elements',
    index_col = 0,
    skiprows = 1)

gradings_data = pd.read_excel(Path("./Input data/") / "test_data_phase_II.xlsx",
    sheet_name='input_rock_gradings',
    index_col = 0,
    skiprows = 2)


# %%
# CALCULATE REQUIRED STONE DIAMETER (Rock armour)

h_list = []
Delta_r_list = []
P_list = []
Sd_allowed_list = []
Hs_list = []
H2_per_list = []
H2_per_Hs_list = []
N_list = []
xi_s_0_2_list = []
Dn50_list = []
h_vdm_shallow_validity_list = []
gamma_beta_stability_list = []
Dn50_oblique_list = []
Dn50_rh_list = []
Dn50_selected_list = []
Delta_c_list = []
Hudson_outcome_list = []
Dn50_concrete_list = []
V_unit_list = []

for Calculation_case in range(1, len(wave_data.index)+1):
    # Calculation_case = 1
    Cross_section_id = wave_data.at[Calculation_case, 'Location']

    # Open project specific parameters
    g           = project_data.at['g'          , 'Value']
    rho_a       = project_data.at['rho_r'      , 'Value']
    rho_w       = project_data.at['rho_w'      , 'Value']
    Cpl_shallow = project_data.at['Cpl_shallow', 'Value']
    Cs_shallow  = project_data.at['Cs_shallow' , 'Value']
    c_beta      = project_data.at['c_beta'     , 'Value']
    c_rh        = project_data.at['c_rh'       , 'Value']
    beta_max    = project_data.at['beta_max'   , 'Value']
    Storm_duration = project_data.at['storm_duration', 'Value']
    Safety         = project_data.at['sf_vdm'        , 'Value']

    #Get info for sea state
    Hm0      = wave_data.at[Calculation_case, 'Hm0']
    wl       = wave_data.at[Calculation_case, 'wl']
    Tm_min_1 = wave_data.at[Calculation_case, 'Tm-1,0']
    dir_wave = wave_data.at[Calculation_case, 'dir_wave']
    Tm_0_2   = wave_data.at[Calculation_case, 'Tm0,2']
    os_bin   = wave_data.at[Calculation_case, 'Offshore bin']
    N        = np.round(Storm_duration*3600/Tm_0_2)


    # Open structure specific parameters
    tana            = cross_section_data.at[Cross_section_id, 'tan_a']
    dir_structure   = cross_section_data.at[Cross_section_id, 'dir_structure']
    z_bed           = cross_section_data.at[Cross_section_id, 'z_bed']
    slope_foreshore = cross_section_data.at[Cross_section_id, 'slope_foreshore']
    P               = cross_section_data.at[Cross_section_id, 'P']
    rh              = cross_section_data.at[Cross_section_id, 'roundhead']


    # Get allowed Sd for the location and calculation case
    LS         = wave_data.at[Calculation_case, 'Limit State']
    Sd_allowed = requirements_data.at['Armour damage limit', LS]

    #Intermediate calculations
    h              = wl-z_bed
    Delta_r          = (rho_a-rho_w)/rho_w
    alpha          = np.arctan(tana)
    waveinfo       = bw.BattjesGroenendijk(Hm0, h, slope_foreshore)
    H2_per         = waveinfo.get_Hp(0.02)
    Hs             = waveinfo.get_Hn(3) #NU BATTJES-GROENENDIJK VOOR Hs UIT Hmo. IS DAT WAT WE WILLEN?

    #OVERRIDE wave characteristics BECAUSE W+B sheet tales Hs = Hm0
    Hs = Hm0
    #H2_per = Hs*1.34
    #user_warning(f"Hs = Hm0, H2_per taken from W+B calculation due to inconsistencies in Battjes Groenendijk")
    xi_s_0_2     = surf_similarity(tana, Hs, Tm_0_2, g)

    # Check validity of Van der Meer shallow
    vdm_shallow_validity = h/Hs

    #Calculate required stone diameter without reduction
    Dn50 = bw.core.vandermeer_shallow(
        Hs, 
        H2_per, 
        Delta_r, 
        P, 
        Sd_allowed, 
        N, 
        xi_s_0_2, 
        alpha, 
        Cpl = Cpl_shallow, 
        Cs = Cs_shallow, 
        safety = Safety
    )
    #Reduce Dn50 with reduction factors from DAR

    #Calculate obliqueness reduction. Function based on SAWP-#3504459-V48-IHS-COA-xxx-CAL_Armour_Stability_under_Waves.XLSM
    beta = calc_beta(dir_structure, dir_wave)
    gamma_beta_stability  = vangent_armour_reduction(beta, c_beta, beta_max)
    Dn50_oblique = Dn50*gamma_beta_stability

    #Dn50_roundhead, no obliqueness correction on roundhead
    Dn50_rh = Dn50*c_rh 

    #Select stone size to apply
    if rh == "Yes":
        Dn50_selected = Dn50_rh
    elif rh == "No":
        Dn50_selected = Dn50_oblique
    else:
        user_warning(f"Incorrect input in Roundhead indication")
        Dn50_selected = 10000

    # Select grading    
    grading = get_class(Dn50_selected, rho_a, gradings_data)
    #It looks like this will only go to the next class if M50>M50_max for a class. Is this the behaviour we want?

    unit = 'Xbloc' #Xbloc or Accropode II

    # Open project specific parameters
    rho_c  = project_data.at['rho_c'      , 'Value'] #For rock armour. Needs adaptation for concrete
    rho_w  = project_data.at['rho_w'      , 'Value']

    # Get info for sea state
    LS     = wave_data.at[Calculation_case, 'Limit State']
    Hm0    = wave_data.at[Calculation_case, 'Hm0']

    # Open cross-section specific parameters
    rh     = cross_section_data.at[Cross_section_id, 'roundhead']

    # Intermediate calculations
    Hs = Hm0 #TO BE CHANGED DEPENDING ON DESIGN APPROACH
    Delta_c = (rho_c-rho_w)/rho_w

    # Create stability case to find outcome for Hudson relationship
    if rh == "Yes":
        stability_case = unit+"-roundhead"
    elif rh == "No":
        stability_case = unit+"-trunk"
    else:
        user_warning(f"Incorrect input in Roundhead identifier")
        stability_case = None

    Hudson_outcome = concrete_element_data.at[LS, stability_case]    

    Dn50_concrete = hudson_fixed_slope(Hs, Hudson_outcome, Delta_c)
    V_unit = Dn50_concrete**3

    h_list.append(h)
    Delta_r_list.append(Delta_r)
    P_list.append(P)
    Sd_allowed_list.append(Sd_allowed)
    Hs_list.append(Hs)
    H2_per_list.append(H2_per)
    H2_per_Hs_list.append(H2_per/Hs)
    N_list.append(N)
    xi_s_0_2_list.append(xi_s_0_2)
    Dn50_list.append(Dn50)
    h_vdm_shallow_validity_list.append(vdm_shallow_validity)
    gamma_beta_stability_list.append(gamma_beta_stability)
    Dn50_oblique_list.append(Dn50_oblique)
    Dn50_rh_list.append(Dn50_rh)
    Dn50_selected_list.append(round(Dn50_selected,4))
    Delta_c_list.append(Delta_c)
    Hudson_outcome_list.append(Hudson_outcome)
    Dn50_concrete_list.append(Dn50_concrete)
    V_unit_list.append(V_unit)

wave_data["h"] = h_list 
wave_data["Delta_r"] = Delta_r_list 
wave_data["P"] = P_list 
wave_data["Sd_allowed"] = Sd_allowed_list 
wave_data["Hs"] = Hs_list 
wave_data["H2_per"] = H2_per_list 
wave_data["H2_per/Hs"] = Hs_list 
wave_data["N"] = N_list 
wave_data["xi_s_0_2"] = xi_s_0_2_list 
wave_data["Dn50"] = Dn50_list 
wave_data["h_vdm_shallow_validity"] = h_vdm_shallow_validity_list 
wave_data["gamma_beta_stability"] = gamma_beta_stability_list 
wave_data["Dn50_oblique"] = Dn50_oblique_list 
wave_data["Dn50_rh"] = Dn50_rh_list 
wave_data["round(Dn50_selected,4)"] = Dn50_selected_list 
wave_data["Delta_c"] = Delta_c_list 
wave_data["Hudson_outcome"] = Hudson_outcome_list 
wave_data["Dn50_concrete"] = Dn50_concrete_list 
wave_data["V_unit"] = V_unit_list 

wave_data.to_excel("wave_data_intermediate_armour_stability.xlsx")

logging.info("Finished intermediate section")

results = []
for location in wave_data.Location.unique():
    location_data = []
    
    location_data.append(location)
    max_LS = max(wave_data[wave_data["Location"] == location]["Dn50"].dropna())
    location_data.append(list(wave_data[wave_data["Dn50"] == max_LS]["Structure"])[0])
    location_data.append(list(wave_data[wave_data["Dn50"] == max_LS]["Limit State"])[0])
    location_data.append(list(wave_data[wave_data["Dn50"] == max_LS]["Offshore bin"])[0])
    location_data.append(list(wave_data[wave_data["Dn50"] == max_LS]["Hm0"])[0])
    location_data.append(list(wave_data[wave_data["Dn50"] == max_LS]["Sd_allowed"])[0])
    location_data.append(max_LS)

    max_LS = max(wave_data[wave_data["Location"] == location]["V_unit"].dropna())
    location_data.append(list(wave_data[wave_data["V_unit"] == max_LS]["Limit State"])[0])
    location_data.append(list(wave_data[wave_data["V_unit"] == max_LS]["Offshore bin"])[0])
    location_data.append(list(wave_data[wave_data["V_unit"] == max_LS]["Hm0"])[0])
    location_data.append(max_LS)

    results.append(location_data)

print(results)
columns = [
    "Location",
    "Structure",
    "Dn50_concrete, LS", 
    "Dn50_concrete, Offshore bin", 
    "Dn50_concrete, Hm0", 
    "Dn50_concrete, Sd_allowed", 
    "Dn50_concrete, max Dn50", 
    "V_unit, LS", 
    "V_unit, Offshore bin",
    "V_unit, Hm0",
    "V_unit, max V_unit",
]
results_df = pd.DataFrame(results, columns=columns)
results_df.to_excel("wave_data_design_armour_stability.xlsx")
logging.info("Finished design section")