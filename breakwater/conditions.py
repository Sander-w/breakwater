import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from .utils.exceptions import InputError, limitstate_warning
from .utils.wave import dispersion
from .core.battjes import BattjesGroenendijk


class LimitState:
    """ Generate LimitState for design

    Define a design limit state, for instance the Ultimate Limit State
    (ULS) or Serviceability Limit State (SLS). The defined limit state[s]
    can then be given as arguments to design one of the supported
    breakwater types, see Table 1 for the required parameters for these
    design classes. A full list of possible arguments can be seen under
    Keyword Arguments.

    Table 1: Required parameters for each design class

    +---------------+------------+------------+------------+------------+
    | Parameter     |    RRM     |  CRM       |  RC        |       CC   |
    +===============+============+============+============+============+
    | h             | x          | x          | x          | x          |
    +---------------+------------+------------+------------+------------+
    | q             | x          | x          | x          | x          |
    +---------------+------------+------------+------------+------------+
    | H13           | o :sup:`1` | o :sup:`1` | o :sup:`1` | o :sup:`1` |
    +---------------+------------+------------+------------+------------+
    | Hm0           | o :sup:`1` | o :sup:`1` | o :sup:`1` | o :sup:`1` |
    +---------------+------------+------------+------------+------------+
    | H2_per        | o :sup:`2` |            |            |            |
    +---------------+------------+------------+------------+------------+
    | Hmax          |            |            | x          | x          |
    +---------------+------------+------------+------------+------------+
    | Tm            | o :sup:`3` |            |            |            |
    +---------------+------------+------------+------------+------------+
    | T13           |            |            | o :sup:`3` | o :sup:`3` |
    +---------------+------------+------------+------------+------------+
    | T_m_min_1     | o :sup:`3` | o :sup:`3` | o :sup:`3` | o :sup:`3` |
    +---------------+------------+------------+------------+------------+
    | Nod           | o :sup:`4` | x          |            |            |
    +---------------+------------+------------+------------+------------+
    | Sd            | o :sup:`4` |            |            |            |
    +---------------+------------+------------+------------+------------+

    | :sup:`1` uses :py:meth:`get_Hs` to get Hs, so at least one of the
      wave heights must be given
    | :sup:`2` if H2_per not specified in the limit state it is computed
      with :py:meth:`get_H2`
    | :sup:`3` wave periods can be transformed in deep water with
      :py:meth:`transform_periods`
    | :sup:`4` either Nod or Sd is required, the other can be computed
      with :py:meth:`Nod` or :py:meth:`Sd`

    .. note::
       The design classes
       :obj:`bw.RockRubbleMound <breakwater.rubble.RockRubbleMound>`,
       :obj:`bw.ConcreteRubbleMound <breakwater.rubble.ConcreteRubbleMound>`
       and :obj:`bw.Caisson <breakwater.caisson.Caisson>` perform
       :sup:`1` and :sup:`2` automatically, but not :sup:`3` and :sup:`4`

    Parameters
    ----------
    h : float
        the water depth [m]
    label : str
        label to identify the LimitState, for instance ULS or SLS

    Keyword Arguments
    -----------------
    q : float
        mean overtopping discharge per meter structure width [l/s per m]
    Hs : float
        the significant wave height [m]
    H13 : float
        mean of the highest 1/3 of the wave heights [m]
    Hm0 : float
        the spectral wave height [m]
    Ho : float
        deep water wave height [m]
    H2_per : float
        the wave height exceeded by 2% of the waves [m]
    Hmax : float
        design wave height for the Goda formula, equal to the mean of
        the highest 1/250 of the wave heights [m].
    H1250 : float
        mean of the highest 1/250 of the wave heights [m], automatically
        interpreted as Hmax.
    Tm : float
        the mean wave period [s]
    T13 : float
         the significant wave period [s]
    Ts : float
         the significant wave period [s]
    T_m_min_1 : float
        :math:`T_{m-1.0}`, the energy wave period [s]
    Tp : float
        the peak period [s]
    Nod : float
        Damage number, used in the formula for the toe stability. Can
        also be computed with :py:meth:`Nod` [-]
    Sd : float
        Damage level parameter, used in Van der Meer formula. Can also
        be computed with :py:meth:`Sd` [-]


    Attributes
    ----------
    conditions : dict
        Dictionary with all parameters
    deep_water : bool
        True if deep water assumption is valid, False if not valid. The
        bool is computed with the Tm, Tp, T_m_min_1 or T13 (in this
        order)
    h : float
        the water depth [m]
    label : str
        identifier of the LimitState
    """
    def __init__(self, h, label, **kwargs):
        """ See help(LimitState) for more info """
        self.conditions = kwargs
        self.h = h
        self.label = label

        if 'H1250' in self.conditions:
            # H1/250 is equal to Hmax thus set Hmax
            self.conditions['Hmax'] = self.conditions['H1250']

        if 'Tm' in self.conditions:
            period = 'Tm'
        elif 'Tp' in self.conditions:
            period = 'Tp'
        elif 'T_m_min_1' in self.conditions:
            period = 'T_m_min_1'
        elif 'Ts' in self.conditions:
            self.conditions['T13'] = self.conditions['Ts']
            period = 'T13'
        elif 'T13' in self.conditions:
            period = 'T13'
        else:
            raise InputError(
                ('No wave period in LimitState, please specify (at least) one'
                 ' of the following wave periods: \'Tm\', \'Tp\', \'Ts\' '
                 '\'T_m_min_1\' or \'T13\'' ))

        if self.h/self.L(period=period) > 0.5:
            self.deep_water = True
        else:
            self.deep_water = False

    def __str__(self):
        return str(self.conditions)

    def __setitem__(self, key, value):
        self.conditions[key] = value

    def __getitem__(self, key):
        return self.conditions[key]

    def check_deep_water(self):
        """ Check if the deep water assumption is valid

        Deep water assumption is valid if h/L > 0.5, where L is computed
        with the dispersion relation using the mean period. Method
        updates the attribute :py:attr:`deep_water` based on the check.

        Raises
        ------
        KeyError
            If Tm is not in the LimitState
        """
        if self.h/self.L(period='Tm') > 0.5:
            self.deep_water = True
        else:
            self.deep_water = False

    def transform_periods(self, scalar):
        """ Transform the wave periods with the Rayleigh characteristics

        Transforms the missing wave periods in the LimitState with the
        characteristics for Rayleigh distributed waves in deep water,
        by using the wave periods which are in the LimitState.

        Table 2: characteristic wave periods (Van den Bos and Verhagen, 2018)

        +-------------------+--------------+
        | Name              |     T/Tp     |
        +===================+==============+
        | :math:`T_{p}`     |      1       |
        +-------------------+--------------+
        | :math:`T_{m}`     | 0.75 to 0.85 |
        +-------------------+--------------+
        | :math:`T_{s}`     | 0.90 to 0.95 |
        +-------------------+--------------+
        | :math:`T_{m-1.0}` |      1.1     |
        +-------------------+--------------+

        Parameters
        ----------
        scalar : float
            Should be between 0 and 1, 0 being the lower bound, 1 the
            upper bound of the characteristics. For instance, scalar of
            1 results in Tm/Tp = 0.85.

        Raises
        ------
        InputError
            If the scalar is not between 0 and 1
        """
        if scalar > 1 or scalar < 0:
            raise ValueError('The scalar should be between 0 and 1')

        if not self.deep_water:
            limitstate_warning(
                'Transforming wave periods is only valid in deep water')

        if 'Tp' not in self.conditions:
            if 'Tm' in self.conditions:
                Tm = self.conditions['Tm']
                self.conditions['Tp'] = Tm/(0.75 + scalar*0.1)
            elif 'T13' in self.conditions:
                T13 = self.conditions['T13']
                self.conditions['Tp'] = T13/(0.9 + scalar*0.05)
            elif 'T_m_min_1' in self.conditions:
                T_m_min_1 = self.conditions['T_m_min_1']
                self.conditions['Tp'] = T_m_min_1/1.1

        if 'Tp' in self.conditions:
            Tp = self.conditions['Tp']
            if 'Tm' not in self.conditions:
                self.conditions['Tm'] = (0.75 + scalar*0.1) * Tp
            if 'T13' not in self.conditions:
                self.conditions['T13'] = (0.9 + scalar*0.05) * Tp
            if 'T_m_min_1' not in self.conditions:
                self.conditions['T_m_min_1'] = 1.1*Tp

    def get_Hs(self, definition):
        """ Get the significant wave height

        Multiple definitions are used for the significant wave height,
        and formulas use different definitions. This method returns the
        desired definition of the significant wave height. For
        overtopping this is Hm0 and H1/3 is used in the Van der Meer
        and Hudson formula.

        .. note::
           If Hm0, Hs or H13 are in the LimitState this method will
           **always** return a significant wave height. The table
           depicts the order in which Hs is returned if a value is
           missing. For example: if the chosen definition is Hm0, then
           first Hm0 is returned, if Hm0 is not included in the
           Limitstate Hs will be retuned, if Hs is also missing H13
           will be returned.

           +------------+-------+-------+-------+
           | definition |  Hm0  |  Hs   |  H13  |
           +============+=======+=======+=======+
           | 1          |  Hm0  |  Hs   |  H13  |
           +------------+-------+-------+-------+
           | 2          |  Hs   |  H13  |  Hs   |
           +------------+-------+-------+-------+
           | 3          |  H13  |  Hm0  |  Hm0  |
           +------------+-------+-------+-------+

        Parameters
        ----------
        definition : {'Hm0', 'Hs', 'H13'}
            definition of the significant wave height to use

        Returns
        -------
        Hs : float
            significant wave height based on the definition

        Raises
        ------
        InputError
            If there is no significant wave height (Hm0, Hs, H1/3) in
            the LimitState
        KeyError
            If the definition is not 'Hm0', 'H13' or 'Hs'
        """
        info = None
        msg = 'No significant wave height in the LimitState'

        if definition == 'Hm0':
            if 'Hm0' in self.conditions:
                H = self.conditions['Hm0']
            elif 'Hs' in self.conditions:
                H = self.conditions['Hs']
                info = f'used Hs instead of {definition} from {self.label}'
            elif 'H13' in self.conditions:
                H = self.conditions['H13']
                info = f'used Hs instead of {definition} from {self.label}'
            else:
                raise InputError(msg)

        elif definition == 'H13':
            if 'H13' in self.conditions:
                H = self.conditions['H13']
            elif 'Hs' in self.conditions:
                H = self.conditions['Hs']
                info = f'used Hs instead of {definition} from {self.label}'
            elif 'Hm0' in self.conditions:
                H = self.conditions['Hm0']
                info = f'used Hm0 instead of {definition} from {self.label}'
            else:
                raise InputError(msg)

        elif definition == 'Hs':
            if 'Hs' in self.conditions:
                H = self.conditions['Hs']
            elif 'H13' in self.conditions:
                H = self.conditions['H13']
                info = f'used H13 instead of {definition} from {self.label}'
            elif 'Hm0' in self.conditions:
                H = self.conditions['Hm0']
                info = f'used Hm0 instead of {definition} from {self.label}'
            else:
                raise InputError(msg)

        else:
            msg = (f'{definition} not a valid key for get_Hs, use \'Hm0\''
                   ', \'H13\' or \'Hs\' instead')
            raise KeyError(msg)

        if info is not None:
            limitstate_warning(info)

        return H

    def L(self, period, deep_water=False):
        """ Compute the wave length with the dispersion relation

        .. math::
           L = L_{o} \\tanh \\left( \\frac{2 \\pi h}{L}\\right)

        in which the deep water wave length is:

        .. math::
           L_{o} = \\frac{g T^{2}}{2 \\pi}

        Parameters
        ----------
        period : str
            keyword argument of a wave period, must be in LimitState
        deep_water : bool, optional, default: False
            if False the dispersion relation will be used, if True
            the deep water wave length will be used. Note that this
            parameter is not equal to the attribute deep_water

        Raises
        ------
        KeyError
            If the specified wave period is not in the LimitState
        """
        T = self.conditions[period]

        if deep_water:
            # user wants deep water wave length
            wave_length = 9.81*T**2/(2*np.pi)
        else:
            # user wants to use the dispersion relation
            wave_length = dispersion(T=T, h=self.h)

        return wave_length

    def s(self, number=None, H=None, T=None):
        """ Compute the fictitious wave steepness

        .. math::
           s_{0} = \\frac{H}{L_{o}}

        the fictitious wave steepness combines the value of the wave
        height at the location of the breakwater with the deep water
        wave length. It is possible to select a number or specify which
        wave height and period must be used for the computation.

        .. note::
           When computing the fictitious wave steepness the significant
           wave height is needed, to get the significant wave height
           form the LimitState the method :py:meth:`get_Hs` is used.

        Parameters
        ----------
        number : {'mean', 'spectral'}, optional, default: None
            definition to use, will automatically select the correct
            definitions for H and T
        H : str, optional, default: None
            wave height in the LimitState
        T : str, optional, default: None
            wave period in the LimitState

        returns
        s_0 : float
            the fictitious wave steepness

        Raises
        ------
        InputError
            If there is no significant wave height (Hm0, Hs, H1/3) in
            the LimitState.
        KeyError
            if the specified wave period is not in the LimitState
        """
        if number == 'mean':
            H = self.get_Hs('H13')
            L = self.L(period='Tm', deep_water=True)
        elif number == 'spectral':
            H = self.get_Hs('Hm0')
            L = self.L(period='T_m_min_1', deep_water=True)
        else:
            H = self.conditions[H]
            L = self.L(period=T, deep_water=True)

        steepness = H/L

        return steepness

    def surf_similarity(
        self, alpha, number=None, H=None, T=None):
        """ Compute the surf similarity parameter

        .. math::
           \\xi = \\frac{tan{\\alpha}}{\\sqrt{s_{o}}}

        Computes the surf similarity parameter, also known as the
        Iribarren number. It is possible to select a number or specify
        which wave height and period must be used for the computation.

        Parameters
        ----------
        alpha : float
            slope of the structure [rad]
        number : {'mean', 'spectral'}, optional, default: None
            definition to use, will automatically select the correct
            definitions for H and T
        H : str, optional, default: None
            wave height in the LimitState
        T : str, optional, default: None
            wave period in the LimitState

        Returns
        -------
        xi : float
            the surf similarity parameter [-]

        Raises
        ------
        InputError
            If there is no significant wave height (Hm0, Hs, H1/3) in
            the LimitState
        KeyError
            If the specified wave period is not in the LimitState
        """
        if H is None and T is None:
            s = self.s(number=number)
            xi = np.tan(alpha)/np.sqrt(s)
        else:
            s = self.s(H=H, T=T)
            xi = np.tan(alpha)/np.sqrt(s)

        return xi

    def get_H2(self, slope_foreshore):
        """ Get the wave height exceeded by 2% of the waves, H2%

        Attempts to return H2% from the defined limit state, if not
        defined in the limit state it is computed with
        :py:class:`BattjesGroenendijk`.

        Parameters
        ----------
        slope_foreshore : float
            the slope of the foreshore [rad]

        Returns
        -------
        H2_per : float
            the wave height exceeded by 2% of the waves [m]

        Raises
        ------
        KeyError
            If Hm0 is not in the LimitState
        """
        if 'H2_per' in self.conditions:
            return self.conditions['H2_per']

        else:
            # get Hm0 and initiate battjes class
            Hm0 = self.conditions['Hm0']
            battjes = BattjesGroenendijk(
                Hm0=Hm0, h=self.h, slope_foreshore=slope_foreshore)

            # add H2% to the LimitState and return the value
            self.conditions['H2_per'] = battjes.get_Hp(0.02)
            return self.conditions['H2_per']

    def Sd(self, G, nv):
        """ Compute the damage level parameter Sd

        This method approximates Sd by using equation 5.150 from the
        Rock Manual (CIRIA, CUR, CETMEF, 2007), and adds Sd to the
        LimitState. The equation to compute Sd is given as:

        .. math::
           S_{d} = \\frac{N_{od}}{G (1 - n_{v})}

        Parameters
        ----------
        G : float
            gradation factor [-]
        nv : float
            porosity of the armour layer [-]

        Returns
        -------
        Sd : float
            damage level parameter [-]

        Raises
        ------
        KeyError
            If Nod is not in the LimitState, use :py:meth:`Nod` instead
        """
        if 'Sd' in self.conditions:
            print('WARNING: This method overwrites the existing Sd with the'
                  ' computed Sd')

        Sd = self.conditions['Nod']/(G*(1-nv))

        self.conditions['Sd'] = Sd

        return Sd

    def Nod(self, G, nv):
        """ Compute the damage number Nod

        This method approximates Nod by using equation 5.150 from the
        Rock Manual (CIRIA, CUR, CETMEF, 2007), and adds Nod to the
        LimitState. The equation to compute Nod is given as:

        .. math::
           N_{od} = G (1 - n_{v}) S_{d}

        Parameters
        ----------
        G : float
            gradation factor [-]
        nv : float
            porosity of the armour layer [-]

        Returns
        -------
        Nod : float
            damage number [-]

        Raises
        ------
        KeyError
            If Sd is not in the LimitState, use :py:meth:`Sd` instead
        """
        if 'Nod' in self.conditions:
            print('WARNING: This method overwrites the existing Nod with the'
                  ' computed Nod')

        Nod = G * (1-nv) * self.conditions['Sd']

        self.conditions['Nod'] = Nod

        return Nod
