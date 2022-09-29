import pandas as pd
import numpy as np

from breakwater.equipment.equipment import Truck, PlateFeeder, Vessel, Excavator, Crane, Barge
from breakwater.utils.breakwater_divide import maxDiff

def install_all_equipments(
        coords,
        equip,
        layer,
        grading,
        mass,
        end_lay,
        area,
        key,
        max_height,
        max_height_xcorner,
        length_top,
        cost_key,
        slope,
        depth_area
    ):
        """
        Function which updates the depth_area dict based on if an equipment can install a section
        Parameters
        ----------
        coords: list
            coordinates of the section
        equip: object
            An equipment object
        layer: str
            Which layer are we installing (e.g. core)
        grading: str
            What is the grading of the layer
        mass: float
            The mass of the material of the layer
        end_lay: float
            Upper y-coordinate of the section
        area: float
            Area of the section
        key: str
            the range of the section (e.g. '0.0-1.0')
        max_height: float
            Maximum height of already build layer
        max_height_xcorner: float
            Corner of heightest build layer
        length_top: float
            Length of highest build layer at max_height
        cost_key: str
            Either CO2 or cost

        Returns
        -------
        tuple
        """
        x, y = list(zip(*coords))
        x, y = np.array(x), np.array(y)
        xcorner_layer = max(coords)[0]
        length_top_sec = maxDiff(x[np.where(y == max(y))])

        # If the equipment is not yet evaluated create a dict, else continue with the already present dict created earlier
        if equip not in depth_area[layer]["Area_yrange"][key]["equipment"].keys():
            depth_area[layer]["Area_yrange"][key]["equipment"][equip] = {}


        # Check to which class the equipment belongs and whether to install
        if isinstance(equip, Truck):
            if equip.install(
                layer=layer,
                grading_layer=grading,
                ymax=max_height,
                section_coords=coords,
                length_top= length_top,
            ):

                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )

                price = area *  equip.get_price(layer=layer, grading_layer=grading, key= cost_key)


                depth_area[layer]["Area_yrange"][key]["equipment"][equip][cost_key] = price
                depth_area[layer]["Area_yrange"][key]["equipment"][equip]['time'] = time

                depth_area[layer]["Area_yrange"][key]["color"] = "g"
                # Is the height of the layer higher than the current max layer? replace top and min_x values
                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, PlateFeeder):
            if equip.install(
                layer=layer,
                grading_layer=grading,
                ymax=max_height,
                section_coords=coords,
            ):

                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )

                price = area *  equip.get_price(layer=layer, grading_layer=grading, key= cost_key)


                depth_area[layer]["Area_yrange"][key]["equipment"][equip][cost_key] = price
                depth_area[layer]["Area_yrange"][key]["equipment"][equip]['time'] = time

                depth_area[layer]["Area_yrange"][key]["color"] = "g"
                # Is the height of the layer higher than the current max layer? replace top and min_x values
                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, Vessel):
            if equip.install(section_coords=coords, layer=layer, grading_layer=grading):
                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )
                price = equip.get_price(layer=layer, grading_layer=grading, key= cost_key) * area


                depth_area[layer]["Area_yrange"][key]["equipment"][equip][cost_key] = price
                depth_area[layer]["Area_yrange"][key]["equipment"][equip]['time'] = time

                depth_area[layer]["Area_yrange"][key]["color"] = "g"

                # Is the height of the layer higher than the current max layer? replace top and min_x values
                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, Excavator):

            if equip.install(
                layer=layer,
                grading_layer=grading,
                ymax=max_height,
                section_coords=coords,
                xmax_top=max_height_xcorner,
                mass=mass,
                length_top=length_top,
                slope= slope,
            ):
                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )
                price = equip.get_price(layer=layer, grading_layer=grading, key= cost_key) * area



                depth_area[layer]["Area_yrange"][key]["equipment"][equip][cost_key] = price
                depth_area[layer]["Area_yrange"][key]["equipment"][equip]['time'] = time

                depth_area[layer]["Area_yrange"][key]["color"] = "g"

                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, Crane):
            if equip.install(
                layer=layer,
                grading_layer=grading,
                ymax=max_height,
                xmax_top=max_height_xcorner,
                length_top= length_top,
                section_coords=coords,
                mass=mass,
                slope = slope
            ):

                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )
                price = equip.get_price(layer=layer, grading_layer=grading, key= cost_key) * area



                depth_area[layer]["Area_yrange"][key]["equipment"][equip][cost_key] = price
                depth_area[layer]["Area_yrange"][key]["equipment"][equip]['time'] = time

                depth_area[layer]["Area_yrange"][key]["color"] = "g"

                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, Barge):

            if equip.install(
                layer=layer,
                grading_layer=grading,
                slope=slope,
                section_coords=coords,
                xmax_top=max_height_xcorner,
                mass=mass,
            ):

                price = (
                    equip.get_price(layer=layer, grading_layer=grading, key= cost_key) * area
                )

                time = area / (
                    equip.get_production_rate(
                        layer=layer, grading_layer=grading
                    )

                )

                depth_area[layer]["Area_yrange"][key]["equipment"][equip][cost_key] = price
                depth_area[layer]["Area_yrange"][key]["equipment"][equip]['time'] = time

                depth_area[layer]["Area_yrange"][key]["color"] = "g"

                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        return max_height, max_height_xcorner, length_top, depth_area

def equipment_dataframe(depth_area):

    d = {
        (i, j): depth_area[i]["Area_yrange"][j]
        for i in depth_area.keys()
        for j in depth_area[i]["Area_yrange"].keys()
    }

    d2 = {}
    area = []
    for layer, section in d.keys():
        n1, n2 = round(float(section.split("-")[0]), 2), round(
            float(section.split("-")[1]), 2
        )
        d2[(layer, f"{n1}-{n2}")] = d[(layer, section)]['equipment']
        area.append(sum(d[(layer, section)]['area']))

    mux = pd.MultiIndex.from_tuples(d2.keys())
    df = pd.DataFrame(list(d2.values()), index=mux)
    df['area'] = area
    cols = [c for c in df.columns if c != 'area']
    cols.insert(0, 'area')
    df = df[cols]
    df = df.applymap(lambda x: None if x == {} else x)

    return df
