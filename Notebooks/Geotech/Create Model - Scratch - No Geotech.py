#%%
# Load packages
from math import nan
from unittest import result
from numpy import NaN
from shapely.geometry.polygon import Polygon
from shapely.geometry import LineString, box

from copy import deepcopy

import datetime

# import geolib
from pathlib import Path

import breakwater as bw

# from geo_dsuite import (
#     AnalysisMethod,
#     DStabilityResults,
#     Headlines,
#     Reflines,
#     SettledStabilityAnalysis,
#     SettledStabilityStage,
#     Materials,
#     SoilLayers,
#     SoilBoundaries,
#     analysis_method,
#     settled_stability,
# )

from typing import List
# from geo_dsuite import soil_layers

# from geo_dsuite.soil_layers import LayerDict

# Define variant to be calculated
variant = "b"

# Breakwater package input
battjes = bw.BattjesGroenendijk(Hm0=4.4, h=15, slope_foreshore=(1, 100))
H2_per = battjes.get_Hp(0.02)
ULS = bw.LimitState(
    h=15,
    Hs=4.5,
    Hm0=4.4,
    H2_per=H2_per,
    Tp=9.4,
    Tm=8.8,
    T_m_min_1=9.7,
    Sd=5,
    Nod=2,
    q=20,
    label="ULS",
)

rho_rock = 2650

NEN = bw.RockGrading(rho=rho_rock)
M50_core = (1 + 1000) / 2
Dn50_core = (M50_core / rho_rock) ** (1 / 3)

RRM = bw.RockRubbleMound(
    slope=(1, 2),
    slope_foreshore=(1, 100),
    rho_w=1025,
    B=5.5,
    N=2100,
    LimitState=ULS,
    Grading=NEN,
    Dn50_core=0.4,
    core_material={"class": "QR 1-1000kg", "Dn50": Dn50_core},
)

variant_polygons = {}
for variant_idx in RRM.variantIDs:
    variant_polygons[variant_idx] = RRM.to_polygon_coordinates(variant_idx)

for variant_idx, variant_dict in variant_polygons.items():
    for key, value in variant_dict.items():
        variant_polygons[variant_idx][key] = Polygon(value)

# Plot variant from Breakwater package
RRM.plot(variant)
