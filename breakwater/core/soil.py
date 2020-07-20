import numpy as np

from ..utils.exceptions import InputError


class Soil:
    """ Define soil

    Define the subsoil on which the breakwater is constructed. This
    class is used in the design classes for making the appropriate 
    geotechnical computations.

    .. note::
       It is currently not possible to use several soil layers in the
       design classes. The soil must therefore be a homogeneous soil.

    Parameters
    ----------
    c : float
        cohesion of the soil [kPa]
    phi : float
        internal friction angle [degrees]
    gamma : float, optional, default: None
        volumetric weight of the soil [kN/m³]
    rho : float, optional, default: None
        the density of the soil [kg/m³]
    n : float, optional, default: None
        porosity of the soil

    Attributes
    ----------
    c : float
        cohesion of the soil [kPa]
    phi : float
        internal friction angle [rad]
    gamma : float
        volumetric weight of the soil [kN/m³]
    gamma_sat : float
        saturated volumetric weight of the soil [kN/m³]
    n : float, optional, default: None
        porosity of the soil
    """

    def __init__(self, c, phi, gamma=None, rho=None, n=None):
        """ See help(Soil) for more info """
        # set input as attributes
        self.c = c
        self.phi = phi * np.pi/180

        # set porosity as attribute
        self.n = n

        # set volumetric weight as attribute
        if gamma is None:
            # check if rho has been given
            if rho is None or n is None:
                raise InputError(
                    ('When the volumetric is not specified, rho and n must '
                     'be specified'))
            else:
                # compute with porosity
                self.gamma = ((1-n)*rho*9.81)/1000
        else:
            # given as input
            self.gamma = gamma

        # set gamma_sat attribute
        self.gamma_sat = None

    def saturated_weight(self, gamma_sat=None, rho_w=None):
        """ Method to add the saturated weight

        Method to add the properties of the soil if it is saturated.
        The saturated soil properties can be added by setting gamma_sat
        or by specifing rho_w. In case of the latter the saturated weight
        is computed with the porosity.

        Parameters
        ----------
        gamma_sat : float, optional, default: None
            saturated volumetric weight of the soil [kN/m³]
        rho_w : float, optional, default: None
            density of water [kg/m³]

        Raises
        ------
        InputError
            if the porosity is not specified and the saturated
            volumetric weight is computed with rho_w
        """
        # check if gamma_sat has been given
        if gamma_sat is not None:
            # set with gamma_sat
            self.gamma_sat = gamma_sat

        else:
            # compute with porosity, check if porosity has been specified
            if self.n is not None:
                # compute saturated volumetric weight
                self.gamma_sat = self.gamma + (self.n*rho_w*9.81)/1000

            else:
                # no porosity has been specified
                raise InputError(
                    ('porosity must be specified when computing the '
                     'saturated volumetric weight'))

    def brinch_hansen(self, p, t, B, L, q, rho_w=None, sat=True):
        """ Brinch-Hansen

        implementation of the Brinch Hansen equation to determine
        the bearing capacity of the soil per unit length. The equation
        is given by (Brinch Hansen, 1970):

        .. math::
           p=i_{c} s_{c} c N_{c}+i_{q} s_{q} q N_{q}+i_{\\gamma}
           s_{\\gamma} \\frac{1}{2} \\gamma B N_{\\gamma}

        In which :math:`i` and :math:`s` are the inclination and shape
        factor, and :math:`N` are dimensionless constants. Note that
        compared to the original equation the depth, base and ground
        inclination factors have been omitted. The latter two because
        the assumption is made that the foundation is never constructed
        at an angle, and the depth factor because the foundation is not
        placed in the soil but on top of the bottom.

        .. note
           the soil is assumed to be homogeneous

        Parameters
        ----------
        t : float
            horizontal stress [kPa]
        p : float
            vertical stress [kPa]
        q : float
            overburden [kPa]
        B : float
            width of the structure [m]
        L : float
            length of the structure [m], set to None if the structure is
            a long structure and the shape factors can be neglected.
        rho_w : float, optional, default: None
            density of water, by default the soil is not submerged [kg/m³]
        sat : bool, optional, default: True
            True if the saturated volumetric weight of the soil must be
            used, False if the dry volumetric weight of the soil must be
            used.

        Returns
        -------
        p : float
            bearing capacity of the foundation [kPa]
        """
        # check density to use
        if sat:
            # check if saturated volumetric weight is given
            if self.gamma_sat is not None:
                # compute effective volumetric weight
                gam = self.gamma_sat - (rho_w*9.81)/1000
            else:
                # not given, but is needed
                raise InputError(
                    ('saturated volumetric weight has not been defined, use '
                     'saturated_weight() to define saturated weight'))
        else:
            # use dry volumetric weight
            gam = self.gamma

        # compute dimensionless constants
        N_q = ((1 + np.sin(self.phi))/(1 - np.sin(self.phi))
                * np.exp(np.pi * np.tan(self.phi)))
        N_c = (N_q - 1)/np.tan(self.phi)
        N_gam = 2 * (N_q - 1) * np.tan(self.phi)

        # compute inclination factors
        i_c = 1 - t / (self.c + p * np.tan(self.phi))
        i_q = i_c**2
        i_gam = i_c**3

        # check if structure is a long structure
        if L is None:
            # shape factors can be neglected, set all shape factors to 1
            s_c, s_q, s_gam = 1, 1, 1
        else:
            # compute shape factors
            s_c = 1 + 0.2 * B/L
            s_q = 1 + np.sin(self.phi) * B/L
            s_gam = 1 - 0.3 * B/L

        # compute capacity of the soil
        p = i_c*s_c*self.c*N_c + i_q*s_q*q*N_q + 0.5*i_gam*s_gam*gam*B*N_gam

        return p
