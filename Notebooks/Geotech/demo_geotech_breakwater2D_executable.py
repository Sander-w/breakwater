#%%
# Load packages
from math import nan
from unittest import result
from numpy import NaN
from shapely.geometry.polygon import Polygon
from shapely.geometry import LineString, box

from copy import deepcopy

import datetime

import geolib
from pathlib import Path

import breakwater as bw

from geo_dsuite import (
    AnalysisMethod,
    DStabilityResults,
    Headlines,
    Reflines,
    SettledStabilityAnalysis,
    SettledStabilityStage,
    Materials,
    SoilLayers,
    SoilBoundaries,
    analysis_method,
    settled_stability,
)

from typing import List
from geo_dsuite import soil_layers

from geo_dsuite.soil_layers import LayerDict

def calculate_stability(breakwater_config, stability_threshold):
    # Plot variant from Breakwater package
    breakwater_config.df.loc[len(breakwater_config.df) - 1].concept.plot("all")

    bw_idx_list = []

    for config_idx in range(0,len(breakwater_config.df)-1):
        # Define config and variant
        # config_idx = len(breakwater_config.df) - 1
        print(f"Configuration no. {config_idx}")
        config = breakwater_config.df.loc[config_idx].concept
        variant = "a"

        # Get coordinates for specific idx from breakwater package input
        variant_polygons = {}
        for variant_idx in config.variantIDs:
            variant_polygons[variant_idx] = config.to_polygon_coordinates(variant_idx)

        # Create polygons from coordinates
        for variant_idx, variant_dict in variant_polygons.items():
            for key, value in variant_dict.items():
                variant_polygons[variant_idx][key] = Polygon(value).buffer(0)

        # Soil parameters input
        materials = Materials.from_excel(
            r"C:\Users\QPQ\Python\Projects\breakwater\Notebooks\Geotech/Materials_with_defaults.xlsx"
        )
        materials.index.name = None
        materials_input = (
            materials.to_soils()
        )

        # Geometry setup
        # Base geometry example
        geom_list = []
        materials_list = []

        geom_list.append(Polygon([(-200, -30), (-200, -20), (200, -20), (200, -30)]))
        materials_list.append(materials_input[9].name)
        geom_list.append(Polygon([(-200, -20), (-200, -10), (200, -10), (200, -20)]))
        materials_list.append(materials_input[5].name)
        geom_list.append(Polygon([(-200, -10), (-200, -5), (200, -5), (200, -10)]))
        materials_list.append(materials_input[4].name)
        geom_list.append(Polygon([(-200, -5), (-200, 0), (200, 0), (200, -5)]))
        materials_list.append(materials_input[2].name)

        # Append breakwater geometry and material to the respective lists
        for variant_idx, variant_dict in variant_polygons.items():
            for key, value in variant_dict.items():
                if variant_idx == variant:
                    geom_list.append(value)
                    materials_list.append(materials.name.loc["Zand_Dijk_Za"])

        # SoilLayers setup
        my_soil_layers = SoilLayers.from_base_input_list(geom_list, materials_list)
        my_soil_layers = my_soil_layers.fix_geometry_consistency()

        # # Plot soillayers
        # my_soil_layers.plot(materials=materials)

        # Headlines setup
        hl1 = LineString([(-50, 10), ((-10 / 2), 10), ((10 / 2 + 2 * 10), 1), (50, 1)])
        my_headlines = Headlines.from_base_input("Freatische lijn (PL1)", hl1, True)
        my_headlines = my_headlines.set_index("name").sort_index()

        # Reflines setup
        rf1 = LineString([(-50, -5), (50, -5)])
        my_reflines = Reflines.from_base_input(
            "Referentie PL1", rf1, "Freatische lijn (PL1)", "Freatische lijn (PL1)"
        )
        my_reflines = my_reflines.set_index("name").sort_index()


        # Analysis method setup
        # Search grid
        search_grid = geolib.models.dstability.analysis.DStabilitySearchGrid(
            bottom_left=geolib.geometry.one.Point(
                label="", id=None, x=-60.0, y=-999.0, z=30, tolerance=0.0001
            ),
            number_of_points_in_x=30,
            number_of_points_in_z=20,
            space=1.0,
        )

        # Analysis method
        my_analysis_method = (
            geolib.models.dstability.analysis.DStabilityBishopBruteForceAnalysisMethod(
                search_grid=search_grid,
                bottom_tangent_line_z=-10,
                number_of_tangent_lines=25,
                space_tangent_lines=0.5,
            )
        )


        # SettledStabilityStage setup
        initial = SettledStabilityStage(
            name="Initial",
            headlines=my_headlines,
            reflines=my_reflines,
            soil_layers=my_soil_layers,
            analysis_method=my_analysis_method,
        )

        # Stability analysis
        analysis = SettledStabilityAnalysis(
            materials=materials,
            stages=[initial],
            path=Path(f"results"),
        )

        my_soil_layers.plot()

        # To settlement model
        m = initial.to_settlement_model(materials=materials)
        m.serialize(
            Path(
                f"results/ScratchSettlementModel_ConfigNr-{config_idx}_Variant-{variant}.sli"
            )
        )

        # To stability model and calculate
        model, mappings = analysis.to_stability_model()
        model.serialize(
            Path(
                f"results/ScratchStabilityModel_ConfigNr-{config_idx}_Variant-{variant}.stix"
            )
        )
        model.execute()

        print(model.get_result(0).FactorOfSafety)

        if model.get_result(0).FactorOfSafety > stability_threshold:
            bw_idx_list.append(config_idx)

    breakwater_config.df = breakwater_config.df[breakwater_config.df.index.isin(bw_idx_list)]

    return breakwater_config




