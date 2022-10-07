import numpy as np
import scipy.special as sc
from scipy.optimize import fsolve


class BattjesGroenendijk:
    """Formulas of Battjes and Groenendijk (2000)

    The wave height distributions on shallow foreshores deviates from
    those in deep water due to the effects of the restricted
    depth-to-height ratio and of wave breaking. Battjes and Groenendijk
    (2000), therefore derived a generalised empirical parameterisations
    based on laboratory data to determine these effects. This allows for
    the computation of all statistical wave heights based on the spectral
    wave height, water depth and slope of the foreshore.

    Parameters
    ----------
    Hm0 : float
        the spectral wave height [m]
    h : float
        the water depth [m]
    slope_foreshore : float or tuple
        the slope of the foreshore [rad], or as tuple formatted as (V,H)

    Attributes
    ----------
    Hrms : float
        the root-mean-square wave height [m]
    Htr_tilde : float
        the transitional wave height, made dimensionless with Hrms [m]
    k1, k2 : floats
        shape parameters of the Composite Weibull Distribution (CWD) [-]
    """

    def __init__(self, Hm0, h, slope_foreshore):
        """ See help(BattjesGroenendijk) for more info """
        # check if slope foreshore is given as a tuple (V, H)
        if isinstance(slope_foreshore, tuple):
            # if so, compute the angle
            slope_foreshore = slope_foreshore[1]/slope_foreshore[0]

        Htr = (0.35 + 5.8*slope_foreshore)*h

        self.Hrms = (0.6725 + 0.2025*(Hm0/h))*Hm0
        self.Htr_tilde = Htr/self.Hrms

        self.k1 = 2
        self.k2 = 3.6

    @staticmethod
    def gammainc_upper(a, x):
        """Upper incomplete gamma function

        Defined as

        .. math::
            \Gamma(a,x) = \int_x^\infty t^{a - 1}e^{-t} dt

        Implemented using scipy.special as gamma(a)*gammaincc(a,x)

        Parameters
        ----------
        a : float
            positive parameter
        x : float
            nonnegative argument

        Returns
        -------
        float
            value of the upper incomplete gamma function
        """
        if a >= 0:
            out = sc.gamma(a)*sc.gammaincc(a, x)
            return out
        else:
            return None

    @staticmethod
    def gammainc_lower(a, x):
        """Lower incomplete gamma function

        Defined as

        .. math::
            \Gamma(a,x) = \int_0^x t^{a - 1}e^{-t} dt

        Implemented using scipy.special as gamma(a)*gammainc(a,x)

        Parameters
        ----------
        a : float
            positive parameter
        x : float
            nonnegative argument

        Returns
        -------
            float
                value of the lower incomplete gamma function
        """
        if a >= 0:
            out = sc.gamma(a)*sc.gammainc(a, x)
            return out
        else:
            return None

    def _solver(self, x):
        """Computes H1~ and H2~ for a given Htr~

        Computes H1~ and H2~ by solving equation 7 from Battjes and
        Groenendijk (2000), and the continuity condition for the
        Composite Weibull Distribution. This continuity conditions is:
        .. math::
            \frac{H_tr}{H_1}^{k_1} = \frac{H_tr}{H_2}^{k_2}

        The values :math:`H_{tr}`, :math:`H_{1}` and :math:`H_{2}` are
        then made nondimensional by dividing them by :math:`H_{rms}`. A
        tilde is added to the names of the values to indicate that they
        are the nondimensional values.

        Parameters
        ----------
        x : list
            initial guess for the values of H1~ and H2~, [float, float]

        Returns
        -------
        list
            containing the continuity equation and equation 7 of
            Battjes and Groenendijk (2000) used to find the values of
            H1~ and H2~ for a given Htr~
        """
        a1 = 2/self.k1 + 1
        a2 = 2/self.k2 + 1

        eq_continuity = (self.Htr_tilde
                         * (self.Htr_tilde/x[0])**(-self.k1/self.k2) - x[1])

        eq7 = (np.sqrt(x[0]**2
                       * self.gammainc_lower(a1,(self.Htr_tilde/x[0])**self.k1)
                       + x[1]**2
                       * self.gammainc_upper(a2,(self.Htr_tilde/x[1])**self.k2))
               - 1)

        return [eq_continuity, eq7]

    def get_Hp(self, P):
        """Computes :math:`H_{P}`

        This method computes the wave height exceeded by P% of the
        waves. It uses scipy.optimize.fsolve to determine H1~ and H2~,
        which are then used to compute :math:`H_{P}`

        Parameters
        ----------
        P : float
            exceedance probability as a fractal [-]

        Returns
        -------
        H_P : float
            wave height exceeded by P\\% of the waves

        Raises
        ------
        ValueError
            If P is entered as a percentage, must be entered as P%/100
        """
        if P >= 1:
            raise ValueError('Do not enter P as a percentage, but P%/100')
        H1_H2 = fsolve(self._solver, [1, 1])

        H_CWD_1 = H1_H2[0]*(-np.log(P))**(1/self.k1)

        if H_CWD_1 <= self.Htr_tilde:
            return H_CWD_1*self.Hrms
        else:
            H_CWD_2 = H1_H2[1]*(-np.log(P))**(1/self.k2)
            return H_CWD_2*self.Hrms

    def get_Hn(self, N):
        """Computes :math:`H_{1/N}`

        This method computes the mean of the highest 1/N-part of the
        wave heights. It uses scipy.optimize.fsolve to determine H1~
        and H2~, which are then used to compute :math:`H_{1/N}`.

        Parameters
        ----------
        N : float
            highest 1/N-part of the wave heights [-]

        Returns
        -------
        H_1/N : float
            mean of the highest 1/N part of the waves

        Raises
        ------
        ValueError
            if N is smaller than 1
        """
        if N >= 1:
            H1_H2 = fsolve(self._solver, [1, 1])

            H1 = H1_H2[0]
            H2 = H1_H2[1]

            a1 = 1/self.k1 + 1
            a2 = 1/self.k2 + 1

            H_N = H1*(np.log(N))**(1/self.k1)

            if self.Htr_tilde >= H_N:
                H_1 = (N*H1*(self.gammainc_upper(a1, np.log(N))
                       - self.gammainc_upper(a1, (self.Htr_tilde/H1)**self.k1))
                       + N*H2
                       * self.gammainc_upper(a2, (self.Htr_tilde/H2)**self.k2))

                return H_1*self.Hrms
            else:
                H_2 = N*H2*self.gammainc_upper(a2, np.log(N))

                return H_2*self.Hrms
        else:
            msg = f'N = {N} has no meaning, please enter a value N >= 1'
            raise ValueError(msg)
