from functools import reduce
import operator
from shapely.geometry import LineString, Point
import numpy as np
import math
import matplotlib.pyplot as plt


def intersect(A, B, C, D):
    try:
        line1 = LineString([Point(A), Point(B)])
        line2 = LineString([Point(C), Point(D)])

        int_pt = line1.intersection(line2)
        return int_pt.x, int_pt.y

    except:
        return None

def GaussianA(x, y):
    A = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
    return A


def clockwise(x, y, slope):
    x = np.array(x)
    xmin = min(x)
    y = np.array(y)
    below_zero = False

    # Transfrom negative values to positive ones
    if any(x[x < 0]):
        x += abs(xmin) + 10
        below_zero = True

    # Get the centre of the coordinates
    coords = list(zip(x, y))
    xc, yc = tuple(
        map(
            operator.truediv,
            reduce(lambda x, y: map(operator.add, x, y), coords),
            [len(coords)] * 2,
        )
    )

    xy_tup = sorted(
        coords,
        key=lambda coord: (
            -135
            - math.degrees(
                math.atan2(*tuple(map(operator.sub, coord, (xc, yc)))[::-1])
            )
        )
        % 360,
    )
    x, y = list(zip(*xy_tup))

    # Transform back if needed
    if below_zero:
        x -= abs(xmin) + 10

    V, H = slope
    x, y = list(x), list(y)

    i = 0
    # It happens that some y-coordinates are on same level but belong to the box below, then remove from the coordinates
    while True:
        if i + 1 == len(x):
            break
        else:
            dy = abs(y[i + 1] - y[i])
            dx = abs(x[i + 1] - x[i])
            # If we move upwards but where not on the slope, remove the coordinates
            if round(dx, 3) > round(dy * (H / V), 3) and dy != 0:
                y.pop(i)
                x.pop(i)
                i = 0
            else:
                i += 1

    return list(x), list(y)

def divide_cross_section(depth_area, slope, plot=False):

        """
        compute the area of a layer for all it's depth ranges

        Parameters
        ----------
        variantID : str
            identifier of the variant, see :py:attr:`variantIDs` for a
            list of all generated variants.

        Returns
        -------
        dict
            coordinates of all layers and the area within certain depth ranges
        """
        # depth_area = self._layers(variantID)
        # slope = self._input_arguments["slope"]
        # self._layers(variantID)
        if plot:
            fig, ax = plt.subplots(figsize=(15, 7.5))
        colorlst = ["g", "b", "y", "r", "orange"]
        p = 0
        for layer, coord in depth_area.items():
            depth_area[layer]["Area_yrange"] = {}

            x_lst = coord["x"]
            y_lst = coord["y"]

            # Create a range list in which we will search for the area's
            n1, n2 = math.ceil(min(y_lst)), math.floor(max(y_lst) + 1)
            range_lst = list(range(n1, n2))
            if len(range_lst) == 1:
                range_lst = [n1, n1]

            xl = coord["x"].copy()
            yl = coord["y"].copy()
            yun, i_un = np.unique(yl, return_index=True)
            xun = np.array(xl)[i_un]
            range_lst.extend(yun)
            xl.extend(xun)
            yl.extend(yun)
            k = 1
            for i in range(len(y_lst) - 1):
                for j in range(len(range_lst)):
                    cor = intersect(
                        (x_lst[i], y_lst[i]),
                        (x_lst[i + 1], y_lst[i + 1]),
                        (min(x_lst), range_lst[j]),
                        (max(x_lst), range_lst[j]),
                    )
                    # if there is an intersection cor is not None
                    if cor != None:
                        x, y = cor

                        xl.insert(i + k, x)
                        yl.insert(i + k, y)
                        k += 1

            # Create a range including the minimum and the maximum
            y_set = sorted(list(set(range_lst)))

            xl, yl = np.array(xl), np.array(yl)
            if plot:
                ax.plot(x_lst, y_lst, "ko")
                ax.fill(x_lst, y_lst, colorlst[p], label=f"{layer}")
                plt.gca().set_aspect("equal", adjustable="box")
                plt.grid()
            p += 1
            for i in range(len(y_set) - 1):
                r1, r2 = y_set[i], y_set[i + 1]
                depth_area[layer]["Area_yrange"][f"{r1}-{r2}"] = {}
                depth_area[layer]["Area_yrange"][f"{r1}-{r2}"]["coordinates"] = []
                # get only the x and y within an y range
                yl2 = yl[(yl >= r1) & (yl <= r2)]
                xl2 = xl[(yl >= r1) & (yl <= r2)]

                depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                    "area"
                ] = []  # Area of a section
                # The layers below can be negative as well as positive with a section interrupted by e.g. the core
                if (
                    layer == "armour"
                    or layer == "underlayer"
                    or layer == "filter layer"
                ):

                    xl2l, xl2r = xl2[xl2 <= 0], xl2[xl2 > 0]
                    yl2l, yl2r = yl2[xl2 <= 0], yl2[xl2 > 0]

                    layers = list(depth_area.keys())
                    layer_below = layers[layers.index(layer) + 1]
                    y_below_max = max(depth_area[layer_below]["y"])

                    # If no y-coordinate is below the top of the layer below -> the section is interrupted in the middle by the later below
                    if not any(yl2 > y_below_max):
                        Al, Ar = 0, 0
                        if len(xl2l) > 0:
                            xl2l, yl2l = clockwise(xl2l, yl2l, slope)
                            if plot:
                                ax.plot(xl2l, yl2l, color="k", linewidth=1)
                            xl2l.append(xl2l[0])
                            yl2l.append(yl2l[0])
                            depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                                "coordinates"
                            ].extend([xl2l, yl2l])
                            Al = GaussianA(xl2l, yl2l)
                            depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                                "area"
                            ].append(Al)
                        if len(xl2r) > 0:
                            xl2r, yl2r = clockwise(xl2r, yl2r, slope)
                            if plot:
                                ax.plot(xl2r, yl2r, color="k", linewidth=1)
                            xl2r.append(xl2r[0])
                            yl2r.append(yl2r[0])
                            depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                                "coordinates"
                            ].extend([xl2r, yl2r])
                            Ar = GaussianA(xl2r, yl2r)
                            depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                                "area"
                            ].append(Ar)

                    else:
                        if len(xl2) > 0:
                            xl2, yl2 = clockwise(xl2, yl2, slope)
                            if plot:
                                ax.plot(xl2, yl2, color="k", linewidth=1)
                            xl2.append(xl2[0])
                            yl2.append(yl2[0])
                            depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                                "coordinates"
                            ].extend([xl2, yl2])
                            A = GaussianA(xl2, yl2)
                            depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                                "area"
                            ].append(A)
                else:
                    if len(xl2) > 0:
                        xl2, yl2 = clockwise(xl2, yl2, slope)
                        if plot:
                            ax.plot(xl2, yl2, color="k", linewidth=1)
                        xl2.append(xl2[0])
                        yl2.append(yl2[0])
                        depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                            "coordinates"
                        ].extend([xl2, yl2])
                        A = GaussianA(xl2, yl2)
                        depth_area[layer]["Area_yrange"][f"{r1}-{r2}"]["area"].append(A)


                depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                    "color"
                ] = "r"  # The starting color of the fill of a section
                depth_area[layer]["Area_yrange"][f"{r1}-{r2}"][
                    "equipment"
                ] = {}  # The used equipment for the section
        if plot:
            plt.legend()

        return depth_area

def maxDiff(a):
        """
        Get the maximum difference between values in an array
        Parameters
        ----------
        a: list

        Returns
        -------
        float
        """

        vmin = a[0]
        dmax = 0
        for i in range(len(a)):
            if a[i] < vmin:
                vmin = a[i]
            elif a[i] - vmin > dmax:
                dmax = a[i] - vmin
        return dmax
