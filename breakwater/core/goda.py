import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from ..utils.exceptions import InputError, user_warning
from ..utils.wave import shoaling_coefficient, dispersion

def goda_wave_heights(h, d, Ho, T, slope_foreshore, Ks=None, factor=1.8):
    """ Compute the design wave height

    Computes the design wave height :math:`H_{1/3}` and
    :math:`H_{max}` with the empirical formulas of Goda (2000)

    Parameters
    ----------
    h : float
        water depth [m]
    d : float
        water depth in front of the caisson, on top of the foundation [m]
    Ho : float
        deep water wave height [m]
    T : float
        wave period, Goda (2000) advises to use :math:`T_{1/3}` [s]
    slope_foreshore : tuple
        slope of the foreshore (V, H). For example a slope of 1:100 is
        defined as (1,100)
    Ks : float, optional, default: None
        non-linear shoaling coefficient. If the shoaling coefficient is
        not specified it is automatically computed with
        :py:func:`shoaling_coefficient` [-]
    factor : float, optional, default: 1.8
        Hmax is a probabilistic quantity. But to avoid possible
        confusion in design, a definite value of Hmax = 1.8 * Ks * H0
        is recommended. The value of 1.8 is, however, a recommendation
        and the user is free to choose another value, such as 1.6 or 2.0
        (Goda, 2000)
    """
    # compute the angle of the foreshore
    slope_foreshore = np.arctan(slope_foreshore[0]/slope_foreshore[1])

    # compute deep water wave length
    L_o = 9.81*T**2/(2*np.pi)

    # check if shoaling coefficient is given
    if Ks is None:
        # if not then compute the non-linear shoaling coefficient
        Ks = shoaling_coefficient(h=d, T=T, H0=Ho)
    else:
        # set the shoaling coefficient to the entered coefficient
        Ks = Ks

    if h/L_o >= 0.2:
        H13 = Ks*Ho
        Hmax = factor*Ks*Ho
    else:
        beta_0 = (0.028*(Ho/L_o)**-0.38
                       * np.exp(20*np.tan(slope_foreshore)**1.5))
        beta_1 = 0.52 * np.exp(4.2*np.tan(slope_foreshore))
        beta_max = max([0.92, (0.32*(Ho/L_o)**-0.29
                               * np.exp(2.4*np.tan(slope_foreshore)))])

        beta_0_star = (0.052*(Ho/L_o)**-0.38
                       * np.exp(20*np.tan(slope_foreshore)**1.5))
        beta_1_star = 0.63 * np.exp(3.8*np.tan(slope_foreshore))
        beta_max_star = max([1.65, (0.53*(Ho/L_o)**-0.29
                                    * np.exp(2.4*np.tan(slope_foreshore)))])

        H13 = min([beta_0*Ho + beta_1*h, beta_max*Ho, Ks*Ho])

        # water depth at a location of 5x H1/3
        hb = h + 5 * np.tan(slope_foreshore) * H13

        Hmax = min([beta_0_star*Ho + beta_1_star*hb,
                    beta_max_star*Ho,
                    factor*Ks*Ho])
    return H13, Hmax

class Goda:
    """ Compute wave pressure with the extended Goda formula (Takahasi, 2002)

    Goda (1992) analysed many of the successful and unsuccessful monolithic
    breakwaters and developed a practical formula that can be used to
    analyse the stability of monolithic breakwaters. The formula
    developed by Goda was not meant to compute the pressures for
    breaking waves (impulsive conditions), therefore Takahasi (2002)
    included an impulsive pressure coefficient in the formula.

    .. warning::
       Goda (1992) advices to avoid impulsive pressures when designing
       monolithic breakwaters.

    Parameters
    ----------
    Hs : float
        mean of the highest 1/3 of the wave heights [m].
    Hmax : float
        design wave height, equal to the mean of the highest 1/250 of
        the wave heights [m].
    h : float
        water depth [m]
    d : float
        water depth in front of the caisson, on top of the foundation [m]
    h_acc : float
        submerged depth of the caisson [m]
    hc : float
        height of the caisson above the water line [m]
    Bm : float
        width of the berm [m]
    T : float
        wave period, Goda (2000) advises to use :math:`T_{1/3}` [s]
    beta : float
        angle between direction of wave approach and a line normal to
        the breakwater [rad]
    rho : float
        density of water [kg/m³]
    slope_foreshore : float
        slope of the foreshore [rad]
    B : float, optional, default: None
        width of the monolithic breakwater [m]
    lambda_ : list, optional, default: [1, 1, 1]
        modification factors of Takahasi (2002) for alternative
        monolithic breakwater. Input must be
        \\lambda_= [:math:`\\lambda_1, \\lambda_2, \\lambda_3`].
    logger : dict, optional, default: None
        dict to log messages, must have keys 'INFO' and 'WARNINGS'

    Attributes
    ----------
    hb : float
        offshore water depth at a distance of five times Hs (=H13) [m]
    L : float
        wave length computed with the dispersion relation [m]
    eta_star : float
        the elevation to which the wave pressure is exerted [m]
    p1, p3, p4 : floats
        representative wave pressure intensities [Pa]
    pu : float
        uplift pressure [Pa]
    h_c_star : float
        elevation to which the wave pressure is exerted on the caisson,
        minimum value of hc and eta_star [m]
    B : float
        width of the monolithic breakwater [m]. None by default so it
        can be computed with :py:meth:`required_width`
    """
    def __init__(
            self, Hs, Hmax, h, d, h_acc, hc, Bm, T, beta, rho,
            slope_foreshore, B=None, lambda_=[1,1,1], logger=None):
        """ See help(Goda) for more info """
        # set dimensions as private variables
        self._h = h
        self._d = d
        self._h_acc = h_acc
        self._hc = hc
        self.rho = rho
        self.B = B

        # water depth at a location of 5x H1/3
        self.hb = h + 5 * np.tan(slope_foreshore) * Hs

        # compute the wave length with the dispersion relation
        self.L = dispersion(T=T, h=self._h)

        # compute wave pressure coefficients (Goda, 2000)
        alpha_1 = (0.6 + 0.5*((4*np.pi*self._h/self.L)
                               / np.sinh(4*np.pi*self._h/self.L))**2)
        alpha_2 = min([(self.hb-self._d)/(3*self.hb)*(Hmax/self._d)**2,
                        2*self._d/Hmax])
        alpha_3 = 1 - self._h_acc/self._h * (1-1/np.cosh(2*np.pi*self._h/self.L))

        # check for impulsive pressures and adjust alpha_2 if needed
        alpha_I = self._impulsive_pressure(Bm=Bm, Hmax=Hmax)
        if alpha_I >= alpha_2:
            alpha_star = alpha_I
            user_warning(
                ('Encountered impulsive conditions in Goda formula, it\'s '
                'recommended to change the dimensions of the structure '
                '(Goda, 2000)'))
        else:
            alpha_star = alpha_2
            if logger is not None:
                logger['INFO'].append(
                    'no impulisve conditions in Goda formula')

        # Compute the elevation to which the wave pressure is exerted
        self.eta_star = 0.75*(1 + np.cos(beta))*lambda_[0]*Hmax

        # Compute the wave pressures
        self.p1 = (0.5*(1 + np.cos(beta))
                   * (lambda_[0]*alpha_1
                      + lambda_[1]*alpha_star*np.cos(beta)**2)
                   * self.rho*9.81*Hmax)
        self.p3 = alpha_3*self.p1
        if self.eta_star > self._hc:
            self.p4 = self.p1*(1-self._hc/self.eta_star)
        else:
            self.p4 = 0
        self.pu = (0.5*(1 + np.cos(beta))
                   * lambda_[2]*alpha_1*alpha_3*self.rho*9.81*Hmax)

        # Determine h_c_star
        self.h_c_star = min([self.eta_star, self._hc])

    @property
    def _width(self):
        """ Private property of the width to check if B can be returned """
        if self.B is None:
            # B was not specified, or has not yet been computed
            # raise error
            raise InputError(
                ('B is not specified as an argument, use required_width() to '
                 'compute the required width of the caisson or specify B as '
                 'an argument'))
        else:
            # B can be returned
            return self.B

    def _impulsive_pressure(self, Bm, Hmax):
        """ Compute impulsive wave pressure coefficient

        Computes the impulsive wave pressure coefficient of the
        extended Goda formula (Takahasi, 2002)
        """
        delta_11 = (0.93*(Bm/self.L-0.12)
                    + 0.36*((self._h-self._d)/self._h-0.6))
        delta_22 = (-0.36*(Bm/self.L-0.12)
                    + 0.93*((self._h-self._d)/self._h-0.6))

        if delta_11 <= 0:
            delta_1 = 20*delta_11
        else:
            delta_1 = 15*delta_11

        if delta_22 <= 0:
            alpha_I1 = np.cos(4.9*delta_22)/np.cosh(delta_1)
        else:
            alpha_I1 = 1/(np.cosh(delta_1)*np.sqrt(np.cosh(3*delta_22)))

        if Hmax <= 2*self._d:
            alpha_I0 = Hmax/self._d
        else:
            alpha_I0 = 2

        alpha_I = alpha_I0*alpha_I1
        return alpha_I

    def _pressure_centroid(self):
        """ Compute the centroid of the pressure

        Reference of the centroid is the lower right of the pressure
        """
        # compute total area
        A_tot = self.p1 * (self._hc + self._h_acc)

        # compute areas to subtract from total area
        A1 = 0.5 * (self.p1 - self.p4) * self._hc
        A2 = 0.5 * (self.p1 - self.p3) * self._h_acc

        # compute centroid of each area, reference is lower right
        y_C_tot = 0.5 * (self._hc + self._h_acc)
        y_C1 = self._h_acc + 2/3 * self._hc
        y_C2 = 1/3 * self._h_acc

        # compute centroid of pressures
        y = (A_tot*y_C_tot - A1*y_C1 - A2*y_C2)/(A_tot - A1 - A2)

        return y

    def _dFv(self, M):
        """ Compute net vertical force """
        # compute net vertical force, down is positive
        return M*9.81 - self.U()

    def effective_width(self, M):
        """ Compute the effective width

        The effective width must be used for geotechnical computations,
        due to the fact that the net vertical force of the caisson is
        eccentric.

        Parameters
        ----------
        M : float
            mass of the caisson [kg]

        Returns
        -------
        float
            the effective width [m]
        """
        return self._width - 2 * self.eccentricity(M)

    def eccentricity(self, M):
        """ Compute the eccentricity of the net vertical force

        Parameters
        ----------
        M : float
            mass of the caisson [kg]

        Returns
        -------
        float
            eccentricity of the net vertical force [m]
        """
        # compute and return eccentricity of the net vertical force
        return self.Ma()/self._dFv(M)

    def P(self):
        """ Compute horizontal force due to the pressure

        Returns
        -------
        float
            horizontal force due to the pressures [Pa]
        """
        P = (0.5*(self.p1+self.p3)*self._h_acc
             + 0.5*(self.p1+self.p4)*self.h_c_star)
        return P

    def Mp(self):
        """ Compute moment at the heel due to the pressure

        Returns
        -------
        float
            moment around the heel due to the horizontal pressures [Nm]
        """
        Mp = (1/6 * (2*self.p1 + self.p3)*self._h_acc**2
              + 0.5*(self.p1 + self.p4)*self._h_acc*self.h_c_star
              + 1/6 * (self.p1 + 2*self.p4)*self.h_c_star**2)
        return Mp

    def U(self):
        """ Compute force due to the uplift pressure

        Returns
        -------
        float
            vertical uplift pressure [Pa]
        """
        return 0.5*self.pu*self._width

    def Mu(self):
        """ Compute moment at the heel due to the uplift

        Returns
        -------
        float
            moment around the heel due to the uplift pressure [Nm]
        """
        return 2/3* self.U() * self._width

    def Ma(self):
        """ Compute the moment around the center of the caisson

        .. warning::
           This method assumes a symmetric caisson

        Returns
        -------
        float
            moment around the center of the caisson [Nm]
        """
        # compute and return moment, clockwise is positive
        return self.U() * self._width/6 + self.P() * self._pressure_centroid()

    def required_mass(
            self, mu, t=0.5, SF_sliding=1.2, SF_turning=1.2, logger=None):
        """ Compute required mass of the monolithic breakwater

        Compute the minimal required mass of the monolithic breakwater
        based on the failure mechanisms sliding and overturning.

        .. math::
           M_{sliding} = \\frac{P SF_{sliding}}{g \\mu} + \\frac{U}{g}
           + \\rho B h'

        .. math::
           M_{turning} = \\frac{M_p SF_{turning}}{g t} + \\frac{M_u}{g t}
           + \\rho B h'

        Parameters
        ----------
        mu : float
            friction factor between the caisson and the foundation [-]
        t : scalar, optional, default: 0.5
            horizontal distance to the centre of gravity [m]
        SF_sliding : float, optional, default: 1.2
            safety factor against sliding. Default value according to
            Goda (2000)
        SF_turning : float, optional, default: 1.2
            safety factor against sliding. Default value according to
            Goda (2000)
        logger : dict, optional, default: None
            dict to log messages, must have keys 'INFO' and 'WARNINGS'

        Returns
        -------
        mass : float
            minimal required mass per meter length to satisfy the
            safety factors [kg/m]
        """
        t = t*self._width

        # Compute forces and moments
        P = self.P()
        Mp = self.Mp()
        U = self.U()
        Mu = self.Mu()

        # Compute the required mass for sliding and turning
        mass_sliding = (P*SF_sliding/(9.81*mu)
                        + U/9.81 + self.rho*self._width*self._h_acc)
        mass_turning = (Mp*SF_turning/(9.81*t)
                        + Mu/(t*9.81) + self.rho*self._width*self._h_acc)

        # normative mass is largest mass of the two
        if mass_sliding >= mass_turning:
            mass_required = mass_sliding
            if logger is not None:
                logger['INFO'].append(
                    ('safety factor for sliding is normative for the '
                     'computation of the mass'))
        else:
            mass_required = mass_turning
            if logger is not None:
                logger['INFO'].append(
                    ('safety factor for overtuning is normative for the '
                     'computation of the mass'))

        return mass_required

    def required_width(
            self, Pc, rho_c, rho_f, rho_w, mu, t=0.5,
            SF_sliding=1.2, SF_turning=1.2, logger=None):
        """ Compute the required width of the monolithic breakwater

        Compute the minimal required width of the monolithic breakwater
        based on the failure mechanisms sliding and overturning.

        Parameters
        ----------
        Pc : float
            contribution of concrete to the total mass of the caisson.
            value between 0 and 1
        rho_c : float
            density of concrete [kg/m³]
        rho_f : float
            density of the fill material, for instance sand [kg/m³]
        rho_w : float
            density of water [kg/m³]
        mu : float
            friction factor between the caisson and the foundation [-]
        t : scalar, optional, default: 0.5
            horizontal distance to the centre of gravity [m]
        SF_sliding : float, optional, default: 1.2
            safety factor against sliding. Default value according to
            Goda (2000)
        SF_turning : float, optional, default: 1.2
            safety factor against sliding. Default value according to
            Goda (2000)
        logger : dict, optional, default: None
            dict to log messages, must have keys 'INFO' and 'WARNINGS'

        Returns
        -------
        B : float
            minimal required width to satisfy the safety factors [m]
        """
        # Compute horizontal force and moment due to this force
        P = self.P()
        Mp = self.Mp()

        # contribution of the fill material to the total mass
        Pf = 1 - Pc

        # compute the mass per meter length and width
        m_acc = (self._hc * (Pf*rho_f + Pc*rho_c)
                 + self._h_acc * (Pf*(rho_f-rho_w) + Pc*(rho_c-rho_w)))

        # compute the minimal required width for sliding and overturning
        B_sliding = 2*SF_sliding*P/(mu*(2*m_acc*9.81 - self.pu))
        B_turning = np.sqrt(SF_turning*Mp
                            / (0.5*m_acc*9.81 - (2*t/3)*self.pu))

        # normative value is largest required width
        if B_sliding >= B_turning:
            B_required = B_sliding
            if logger is not None:
                logger['INFO'].append(
                    ('safety factor for sliding is normative for the '
                     'computation of the width'))
        else:
            B_required = B_turning
            if logger is not None:
                logger['INFO'].append(
                    ('safety factor for overtuning is normative for the '
                     'computation of the width'))

        # update attribute, first check if one wass specified
        if self.B is not None:
            # B was not None, show warning
            user_warning(
                f'B was given as {self.B}, now changed to {B_required}')

        # update B
        self.B = B_required

        return B_required

    def bearing_pressure_width(self, B1, Pc, rho_c, rho_fill, pe_max, t=0.5):
        """ Compute the required width for the bearing pressure

        Method uses fsolve from scipy to compute the width that
        satisfy the maximum bearing pressure of the foundation.

        Parameters
        ----------
        B1 : float
            first estimate for the width [m]
        Pc : float
            contribution of concrete to the total mass of the caisson.
            value between 0 and 1
        rho_c : float
            density of concrete [kg/m³]
        rho_f : float
            density of the fill material, for instance sand [kg/m³]
        pe_max : float
            maximum value of the bearing pressure at the heel of the
            caisson. Goda (2000) advises a value between 400 and 500
            kPa.
        t : scalar, optional, default: 0.5
            horizontal distance to the centre of gravity [m]

        Returns
        -------
        float
            required width for the maximum bearing pressure [m]
        """
        # define lambda function for fsolve and compute new width
        func = lambda B: (self.bearing_pressure(
                B=B, Pc=Pc, rho_c=rho_c, rho_fill=rho_fill)
            - pe_max*1000)
        B = fsolve(func, B1, xtol=0.02)[0]

        # update attribute
        self.B = B

        return B

    def bearing_pressure(self, Pc, rho_c, rho_fill, t=0.5, B=None):
        """ compute the bearing pressure at the heel

        Method to compute the bearing pressure at the heel of the
        caisson.

        .. math::
           p_{e}=\\left\\{\\begin{array}{ll} \\frac{2 W_{e}}{3 t_{e}} & :
           t_{e} \\leq \\frac{1}{3} B \\\ \\frac{2 W_{e}}{B}\\left(2-3
           \\frac{t_{e}}{B}\\right) & : t_{e}>\\frac{1}{3} B
           \\end{array}\\right.

        in which:

        .. math::
           t_{e}=\\frac{M_{e}}{W_{e}}, \\quad M_{e}=M g t-M_{U}-M_{p},
           \\quad W_{e}=M g-U

        Parameters
        -------
        Pc : float
            contribution of concrete to the total mass of the caisson.
            value between 0 and 1
        rho_c : float
            density of concrete [kg/m³]
        rho_f : float
            density of the fill material, for instance sand [kg/m³]
        t : scalar, optional, default: 0.5
            horizontal distance to the centre of gravity [m]
        B : float, optional, default: None
            width of the monolithic breakwater [m], used with
            :py:meth:`bearing_pressure_width` to compute the required
            width to satisfy the maximum bearing pressure.

        Returns
        -------
        pe : float
            the bearing pressure at the heel of the caisson [Pa]
        """
        if B is None:
            # get the width
            B = self._width

        # compute the mass of the structure
        M = self.mass(Pc=Pc, rho_c=rho_c, rho_fill=rho_fill)

        W = (M - self.rho*B*self._h_acc) * 9.81

        We = W - self.U()
        Me = W*t*B - self.Mu() - self.Mp()
        te = Me/We

        # compute bearing pressure
        if te <= B/3:
            pe = 2*We/(3*te)
        else:
            pe = 2*We/B * (2 - 3*te/B)

        return pe

    def mass(self, Pc, rho_c, rho_fill):
        """ Compute the mass of the caisson

        Parameters
        ----------
        Pc : float
            contribution of concrete to the total mass of the caisson.
            value between 0 and 1
        rho_c : float
            density of concrete [kg/m³]
        rho_f : float
            density of the fill material, for instance sand [kg/m³]

        Returns
        -------
        m : float
            minimal required mass per meter length to satisfy the
            safety factors [kg/m]
        """
        # compute area
        A = self._width * (self._h_acc + self._hc)

        # compute the mass
        M = A*Pc*rho_c + A*(1-Pc)*rho_fill

        return M

    def plot(self):
        """ Plot pressure distribution

        Plots the pressure distribution together with the dimensions of
        the monolithic breakwater.

        .. warning::
           Do not read the dimensions of the monolithic breakwater from
           the axes of the figure. The correct dimensions of the
           monolithic breakwater can be read from the figure.
        """
        # get the width
        B = self._width

        p = np.array([self.p3,
                      self.p1,
                      self.p1-self.p1*self._hc/self.eta_star])/1000
        y = np.array([0, self._h_acc, self._h_acc+self._hc])
        x = np.array([0, B])
        pu = np.array([self.pu, 0])/1000

        # scale the size of the caisson with the pressure
        scale = np.max(p)/np.max(y)
        y = y*scale
        x = x*scale

        # plot the caisson and wlev
        plt.vlines(x=0, ymin=0, ymax=-y[2], linewidth=2)
        plt.vlines(x=-x[1], ymin=0, ymax=-y[2], linewidth=2)
        plt.hlines(y=0, xmin=0, xmax=-x[1], linewidth=2)
        plt.hlines(y=-y[2], xmin=0, xmax=-x[1], linewidth=2)

        plt.hlines(y=-self._h_acc*scale, xmin=np.max(p)*1.3, xmax=0, color='b')

        # plot dimensions of the caisson
        plt.text(
            x=(-x[1]/2)*1.1, y=(-y[2]/2)*1.03, s=f'B = {np.round(B, 1)} m')
        plt.arrow(x=-x[1]/2, y=-y[2]/2, dx=x[1]/2, dy=0, head_width=4,
                  head_length=4, color='gray', length_includes_head=True)
        plt.arrow(x=-x[1]/2, y=-y[2]/2, dx=-x[1]/2, dy=0, head_width=4,
                  head_length=4, color='gray', length_includes_head=True)

        plt.text(x=(-x[1]/2)*0.9, y=(-y[2]/2)*1.1,
                 s=f'L = {np.round(y[2]/scale, 2)} m', rotation=90)
        plt.arrow(x=-x[1]/2, y=-y[2]/2, dx=0, dy=-y[2]/2, head_width=4,
                  head_length=4, color='gray', length_includes_head=True)
        plt.arrow(x=-x[1]/2, y=-y[2]/2, dx=0, dy=y[2]/2, head_width=4,
                  head_length=4, color='gray', length_includes_head=True)

        # plot pressure distributions
        plt.hlines(y=0, xmin=0, xmax=p[0], color='deepskyblue')
        plt.hlines(y=-y[2], xmin=0, xmax=p[2], color='deepskyblue')
        plt.plot(p, -y, color='deepskyblue')

        plt.vlines(x=0, ymin=0, ymax=pu[0], color='deepskyblue')
        plt.plot(-x, pu, color='deepskyblue')

        # hide negative ticks
        xticks = [xtick for xtick in plt.gca().get_xticks() if xtick >= 0]
        yticks = [ytick for ytick in plt.gca().get_yticks() if ytick >= 0]

        plt.gca().set_xticks(xticks)
        plt.gca().set_yticks(yticks)

        # invert axes
        plt.xlim(np.max(p)*1.3, -np.max(x)-5)
        plt.ylim(np.max(pu)+5, -np.max(y)-5)

        # add title and label
        plt.title('Pressure distributions computed with Goda (2000)')
        plt.xlabel('Pressure [kPa]')
        plt.ylabel('Pressure [kPa]')

        # add grid and show plot
        plt.grid()
        plt.show()
