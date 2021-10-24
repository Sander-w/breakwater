from functools import reduce
import operator
from shapely.geometry import LineString, Point
import numpy as np
import math


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
