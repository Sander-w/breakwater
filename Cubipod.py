import breakwater as bw


class Cubipod(bw.ConcreteArmour):
    """Custom Armour Unit class"""

    # define the init method, note that it is not required to give a
    # keyword argument.
    def __init__(self, rho=2400):
        """See help(CustomArmourUnit) for more info"""
        # define the armour units in the init method
        units = {
            1.0: {"D": 1.0, "h": 1.0, "Vc": 0.59},
            1.5: {"D": 1.15, "h": 1.1, "Vc": 0.68},
            3.0: {"D": 1.43, "h": 1.4, "Vc": 0.85},
            5.1: {"D": 1.72, "h": 1.7, "Vc": 1.02},
            8.1: {"D": 2.01, "h": 2.0, "Vc": 1.18},
            12.1: {"D": 2.29, "h": 2.3, "Vc": 1.35},
            23.6: {"D": 2.87, "h": 2.9, "Vc": 1.69},
            46.1: {"D": 3.59, "h": 3.6, "Vc": 2.12},
            79.7: {"D": 4.30, "h": 4.3, "Vc": 2.54},
        }

        # call the init method of bw.ConcreteArmour with super and pass
        # the following arguments: kd factor for Hudson formula, name of
        # the armour unit for overtopping and rho for the density
        super().__init__(kd=10, units=units, name="Cubipods", rho=rho)

    # this is an optional step. It is possible to define a method
    # to compute a correction factor to be applied on the computed
    # required nominal diameter. Method must allow for kwargs input
    # as several parameters are passed to the method, it is not
    # necessary to use all of the passed parameters
    def correction_factor(self, h, Hs, Rc, **kwargs):
        # it is then possible to build logic in order to determine the
        # correction factor, for example:
        # define list to store more correction factors
        correction = []

        # check for water depth wave height relation
        if h > 2.5 * Hs:
            correction.append(1.3)
        else:
            pass

        # check for low crested structure
        if Rc / Hs < 1:
            correction.append(1.5)
        else:
            pass

        # check if any correction factors have been added
        if any(correction):
            # return maximum value
            return max(correction)
        else:
            # return 1, in other words, no correction factor
            # make sure a correction_factor is always returned
            return 1
