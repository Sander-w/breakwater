import numpy as np
import matplotlib.pyplot as plt

from ..utils.exceptions import user_warning, InputError


class Bishop:
    """ Bishop slip circles

    Class for the computation of the factor of safety against slip
    failure. This factor of safety is computed with the following
    equation (Verruijt, 2012):

    .. math::
       F=\\frac{\\sum \\frac{c+(\\gamma h-p) \\tan \\phi}{\\cos \\alpha(
       1+\\tan \\alpha \\tan \\phi / F)}}{\\sum \\gamma h \\sin \\alpha}

    The top of the soil is defined by one point, point2, see the Figure.
    All slip circles must go through (0,0) and point2, this means that
    if a custom slip circle is defined this circle must go through both
    points. If a custom slip circle is not given, slip circles are
    automatically generated. These circles are generated with an interval
    of :py:obj:`y_step`.

    .. figure:: _figures/bishop-def.png
       :align: center
       :alt: definition of the input parameters for the Bishop class

    Parameters
    ----------
    point2 : tuple
        x and y coordinate of point2
    y_step : float, optional, default: 1
        step size of the y coordinate for generating circles, used if a
        custom circle is not defined
    wlev : float, optional, defautl: None
        y coordinate of the water level, by default the water level is
        None, which means that there is no water
    SlipCircle : :py:class:`SlipCircle`
        user defined custom slip circle, circle must be defined with
        :py:class:`SlipCircle`

    Attributes
    ----------
    x2, y2 : float
        x and y coordinate of point2
    wlev : float
        y coordinate of the water level
    circles : dict
        dictionary with the centre and radius of the circle. Slices are
        added to the dict when generated during the computation.
    layers : dict
        dictionary with the defined layers
    normative : int
        id of the normative :py:class:`SlipCircle`
    """

    def __init__(self, point2, y_step=1, wlev=None, SlipCircle=None):
        """ See help(Bishop) for more info """
        # set point2 and y_step as attributes
        self.x2 = point2[0]
        self.y2 = point2[1]
        self.wlev = wlev

        # check if a custom circle has been specified
        if SlipCircle is None:
            # no circle given, generate all circles
            self.circles = self._make_circles(y_step)
        else:
            # set custom circle as attribute with id 0
            self.circles = {1: SlipCircle}

        # make attribute for the layers
        self.layers = {}

        # set attribute to store id of normative slip circle
        self.normative = None

    def _make_circles(self, y_step):
        """ Compute all possible circles

        This function computes all possible circles going through two
        points, the first is (0,0) and the second must be specified as
        an argument. The possible circles are computed with the following
        constraints:

        .. math::
           y \\geq y_{2}
           0 \\leq x \\leq x_{2}

        Parameters
        ----------
        y_step : float, optional, default: 1
            step size with which the y coordinate is increased

        Returns
        -------
        dict
            dictionary with the centre coordinate, radius and angle of the
            arc for all possible circles
        """
        # make dict to store generated circles
        circles = {}

        # set start value for y of the centre of the circle
        y = self.y2

        # set variable for the circle id
        id = 1

        # generate circles in the domain
        while True:
            # compute x coordinate of the centre of the circle
            x = (self.x2**2 + self.y2**2 - 2*self.y2*y)/(2*self.x2)

            # compute radius of the circle
            r = np.sqrt(x**2+y**2)

            # compute angle of the arc
            chord = np.sqrt(self.x2**2 + self.y2**2)

            # make circle object and add to dict
            circles[id] = SlipCircle(centre=(x, y), r=r)

            # check if x is negative, break criteria
            if x <= 0:
                # break loop
                break
            else:
                # increase y with the step size and increment i
                y += y_step
                id += 1

        return circles

    def add_layer(
            self, gamma, c, phi, name, ymin=None, ymax=None, gamma_sat=None):
        """ Add layer to the soil

        Add a soil layer for the computation. In case of a homogeneous
        soil the arguments :py:obj:`ymin` and :py:obj:`ymax` do not have
        to be specified. However, if the soil consists of several layers
        these parameters must be given for each layer.

        .. warning::
           soil layers must be added sequentially from the lowest layer
           to the highest layer.

        Parameters
        ----------
        gamma : float
            volumetric weight of the material [kN/m³]
        c : float
            cohesion of the material [kPa]
        phi : float
            internal friction angle [deg]
        name : str
            name of the layer
        ymin : float, optional, default: None
            y coordinate of the start of the layer, required if more
            than 1 layer is added
        ymax : float, optional, default: None
            y coordinate of the end of the layer, required if more
            than 1 layer is added
        gamma_sat : float, optional, default: None
            saturated volumetric weight of the material [kN/m³]. Only
            required if a :py:attr:`wlev` is specified.

        Raises
        ------
        InputError
            if more than 1 layer is specified and ymin and/or ymax has
            not been specified or gamma_sat is not given when a wlev has
            been specified
        """
        # convert phi to radians
        phi = phi*np.pi/180

        # check if name is already in
        if name in self.layers.keys():
            # show warning that layer is overwritten
            user_warning(
                (f'{name} already a layer, layer has been overwritten by the'
                  ' specified layer'))

        # check if a wlev has been defined, and a gamma_sat
        if self.wlev is not None and gamma_sat is None:
            raise InputError(
                'gamma_sat must be given when a wlev is specified')

        # add to layers attribute
        self.layers[name] = {
            'gamma': gamma, 'gamma_sat': gamma_sat, 'c': c, 'phi': phi,
            'y_coord': (ymin, ymax)}

        # check if layers are correctly added
        if self.layers:
            # a layer has already been added, iterate over the layers
            # since layers must be added sequentially
            for i, (layer, params) in enumerate(self.layers.items()):
                # get y coordinates of the layer
                ymin_current, ymax_current = params['y_coord']

                # start checking for more than 1 layer
                if i > 0:
                    # next layer, check if ymax of the previous layer
                    # is equal to ymin of the this layer and ymax is
                    # larger since layer must be on top of the previous
                    # layer
                    if (ymax_previous == ymin_current
                            and ymax_current > ymax_previous):
                        # correct
                        pass
                    else:
                        raise InputError(
                            (f'layer {layer} has invalid coordinates, layers '
                              'must be added sequentially'))

                # save coordinates of the current layer
                ymin_previous, ymax_previous = ymin_current, ymax_current

    def compute(self, num_slices, max_iter=50, ftol=0.05, gamma_w=None):
        """ Compute the factor of safety

        Method to compute the factor of safety for all generated
        circles, or the specified circle. Updates :py:attr:`normative`
        with the id of the normative slip circle, i.e. the slip circle
        with the lowest factor of safety.

        Parameters
        ----------
        num_slices : int
            number of slices
        max_iter : int, optional, default: 50
            maximum number of iterations
        ftol : float, optional, default: 0.05
            break criterium, when the change between the previous and
            the current factor of safety is below the change the
            iteration is ended
        gamma_w : float, optional, default: None
            volumetric weight of water [kN/m³]
        """
        # check if max_iter is an int
        if not isinstance(max_iter, int):
            raise InputError('max_iter must be an int')

        # set variable to track normative F
        F_norm = 10*10

        # iterate over the circles
        for id, SlipCircle in self.circles.items():
            # check if slices have alreay been made
            if SlipCircle.slices is None:
                # slices not have not yet been made, so make slices
                SlipCircle.make_slices(
                    num_slices=num_slices, point2=(self.x2, self.y2),
                    layers=self.layers)

            # set first estimate for F and iteration tracker
            F = 1
            i = 0

            # iteratively compute the factor of safety
            while True:
                # set strength and load variable
                strength, load = 0, 0

                # iterate over the slices
                for slice, coords in SlipCircle.slices.items():
                    # check slice number
                    if slice == 0:
                        # only include pressure for the first slice
                        pressure = True
                    else:
                        # other slides no pressure, otherwise doubling
                        pressure = False

                    # compute load and strength
                    load += self._load(coords['alpha_s'], coords['h'])
                    strength += self._strength(
                        alpha_s=coords['alpha_s'], heights=coords['h'], F=F,
                        pressure=pressure, gamma_w=gamma_w)

                # # compute factor of safety
                F_new = strength/load

                # check for break criteria
                if np.abs(F_new - F) <= ftol or i == (max_iter-1):
                    # set F as attribute of the slip circle
                    SlipCircle.F = F_new

                    # check if computed F is smaller than current normative
                    if F_new <= F_norm:
                        # update F_norm and set id as normative
                        F_norm = F_new
                        self.normative = id

                    # break computation
                    break

                # update F and increment i
                F = F_new
                i += 1

    def _load(self, alpha_s, heights):
        """ Compute the load on the slice

        Method computes the loading force of a slice with the following
        equation:

        .. math::
           \\text{load} = \\gamma h \\sin \\alpha

        Parameters
        ----------
        alpha_s : float
            slip angle of the slice [rad]
        heights : dict
            dictionary with the height of each layer

        Returns
        -------
        float
            load of the slice
        """
        # set variable for adding weights
        W = 0

        # iterate over the layers
        for layer, h in heights.items():
            # check if there is a wlev
            if self.wlev is not None:
                # check if layer is above or below the wlev
                if self.wlev >= np.max(h):
                    # entire layer is below the water level
                    W += self.layers[layer]['gamma_sat']*(np.max(h)-np.min(h))

                elif np.min(h) >= self.wlev:
                    # entire layer is above the wlev
                    W += self.layers[layer]['gamma']*(np.max(h)-np.min(h))

                else:
                    # layer is partly in the water, compute wet and dry h
                    h_dry = np.max(h) - self.wlev
                    h_wet = (np.max(h)-np.min(h)) - h_dry

                    # compute weight of the layer
                    W += (self.layers[layer]['gamma_sat']*h_wet
                          + self.layers[layer]['gamma']*h_dry)

            else:
                # not wlev specified, compute the load
                W += self.layers[layer]['gamma']*(np.max(h)-np.min(h))

        return W*np.sin(alpha_s)

    def _strength(self, alpha_s, heights, F, pressure, gamma_w=None):
        """ Compute the strength of a slice

        Method computes the loading force of a slice with the following
        equation:

        .. math::
           \\text{strength} = \\sum \\frac{c+(\\gamma h-p) \\tan \\phi}
           {\\cos \\alpha(1+\\tan \\alpha \\tan \\phi / F)}

        Parameters
        ----------
        alpha_s : float
            slip angle of the slice [rad]
        heights : dict
            dictionary with the height of each layer
        F : float
            factor of safety [-]
        pressure : bool
            True if the pressure must be included in the computation,
            False if not
        gamma_w : float, optional, default: None
            volumetric weight of water [kN/m³]

        Returns
        -------
        float
            strength of the slice
        """
        # set variable for adding weights
        W = 0

        # check if a wlev is defined
        if self.wlev is not None and pressure:
            # check if gamma_w has been given
            if gamma_w is None:
                # no volumetric weight of water, raise error
                raise InputError(
                    'gamma_w must be specified when computing with a wlev')

            # get top and bottom coordinate of the slice
            ytop = np.max(list(heights.values()))
            ybottom = np.min(list(heights.values()))

            # check if top is above the wlev
            if ytop >= self.wlev:
                # compute pressure in the middle of the slice
                p = 0.5*(self.wlev-ybottom)*gamma_w
            else:
                # top of the slice is not above the wlev
                # compute pressure in the middle of the slice
                p = gamma_w*(self.wlev-ytop+0.5*(ytop-ybottom))

        else:
            # no wlev, thus pressure is zero
            p = 0

        # iterate over the layers
        for i, (layer, h) in enumerate(heights.items()):
            # check if first layer
            if i == 0:
                # first layer of the slice offers the resistance
                # against sliding, thus save parameters of this layer
                params = self.layers[layer]

            # check if there is a wlev
            if self.wlev is not None:
                # check if layer is above or below the wlev
                if self.wlev >= np.max(h):
                    # entire layer is below the water level
                    W += self.layers[layer]['gamma_sat']*(np.max(h)-np.min(h))

                elif np.min(h) >= self.wlev:
                    # entire layer is above the wlev
                    W += self.layers[layer]['gamma']*(np.max(h)-np.min(h))

                else:
                    # layer is partly in the water, compute wet and dry h
                    h_dry = np.max(h) - self.wlev
                    h_wet = (np.max(h)-np.min(h)) - h_dry

                    # compute weight of the layer
                    W += (self.layers[layer]['gamma_sat']*h_wet
                          + self.layers[layer]['gamma']*h_dry)

            else:
                # not wlev specified, compute the weight of the layer
                W += params['gamma']*(np.max(h)-np.min(h))

        # compute the strength
        A = params['c'] + (W - p)*np.tan(params['phi'])
        B = np.tan(alpha_s)*np.tan(params['phi'])/F
        C = np.cos(alpha_s)*(1 + B)
        strength = A/C

        return strength

    def plot(self, id=None, show_slices=False):
        """ Plot slip circle(s)

        Method to plot the slip circle(s)

        Parameters
        ----------
        id : int, optional, default: None
            id of the slip circle to plot, by default all slip circles
            are plotted
        show_slices : bool, optional, default: False
            if the slices must be shown, only available after the
            computation as the slices are not generated before the
            computation. The default value is False, meaning that the
            slices will not be plotted.
        """
        # create figure
        fig, ax = plt.subplots()

        # define coordinates of the boundary
        x = [-self.x2, 0, self.x2, 2*self.x2]
        y = [0, 0, self.y2, self.y2]

        # plot the boundary
        ax.plot(x,y, zorder=20, color='k', lw=2)

        # check if custom id has been specified
        if id is not None:
            # plot one circle
            ax = self.circles[id]._make_figure(ax, show_slices)

            # compute ymax and ymin
            ymax = self.circles[id].xy[1] + self.circles[id].r
            ymin = self.circles[id].xy[1] - self.circles[id].r

            # check if ymin is zero
            if ymin >= -1 or ymin <= 1:
                ymin = -5
                
            # set title
            title = f'Slip circle id={id}'

        else:
            # set variables for ymax and ymin
            ymax, ymin = 0, 0

            # iterate over the generated circles to plot all circles
            for id, SlipCircle in self.circles.items():
                # plot circle
                ax = SlipCircle._make_figure(ax, show_slices=False)

                # compute ymax and ymin
                ymax1 = self.circles[id].xy[1] + self.circles[id].r
                ymin1 = self.circles[id].xy[1] - self.circles[id].r

                # check if largest value is larger than current largest
                if ymax1 > ymax:
                    ymax = ymax1

                # check if smallest value is smaller than current smallest
                if ymin1 < ymin:
                    ymin = ymin1

            # set title
            title = 'Slip circles'

        # check if more than one layer has been defined
        if len(self.layers) > 1:
            # plot the layers
            colors = ['g', 'r', 'c', 'm', 'y']
            for i, (layer, params) in enumerate(self.layers.items()):
                # check if out of range for colors
                if i == len(colors):
                    # reset i
                    i = 0

                # plot layer
                ax.axhline(
                    y=params['y_coord'][1], label=f'top of {layer}', ls='--',
                    color=colors[i])

            plt.legend()

        # check if a wlev has been specified
        if self.wlev is not None:
            # add wlev
            plt.axhline(y=self.wlev, color='blue')

        # set layout
        # ax.set_xlim(np.min(x), np.max(x))
        ax.set_ylim(ymin*1.1, ymax*1.1)
        ax.set_xlabel('x coordinate [m]')
        ax.set_ylabel('y coordinate [m]')
        ax.set_title(title)
        ax.grid()
        fig.gca().set_aspect('equal', adjustable='box')
        fig.tight_layout()
        plt.show()


class SlipCircle:
    """ Define a slip circle

    Define a custom slip circle to use in :py:class:`Bishop`

    Parameters
    ----------
    centre : tuple
        x and y coordinate of the centre of the circle
    r : float
        radius of the circle

    Attributes
    ----------
    xy : tuple
        x and y coordinate of the centre of the circle
    r : float
        radius of the circle
    F : float
        computed factor of safety
    slices : dict
        dictionary with the coordinates, height and slip angle of the
        slices
    """

    def __init__(self, centre, r):
        """ See help(SlipCircle) for more info """
        # set input as attributes
        self.xy = centre
        self.r = r

        # set attribute for factor of safety and slices
        self.F = None
        self.slices = None

    def __repr__(self):
        return f'SlipCircle(centre={self.xy}, r={np.round(self.r, 2)})'

    def __str__(self):
        return (f'Slip cricle with centre {self.xy} and a radius of '
                f'{np.round(self.r, 2)}')

    @staticmethod
    def _compute_heights(y1, y2, y3, y4, layers):
        """ Compute the height of the slice

        Method to compute the height of a slice. Note that the following
        condition must hold in the input:

        .. math::
        y_4 > y_3 > y_2 > y_1

        Parameters
        ----------
        y1, y2 : float
        lower y coordinates of the slice, where y2 must be larger
        than y1
        y3, y4 : float
        higher y coordinates of the slice, where y4 must be larger
        than y3

        """
        # create dict to store the heights
        height = {}

        # compute y coordinates of the middle of the circle
        ym_low = y1 + 0.5*(y2-y1)
        ym_top = y3 + 0.5*(y4-y3)

        # check if more than one layers has been defined
        if len(layers) > 1:
            # set variable to track progress through the layers
            # starting value of this variable is ym_low
            y = ym_low

            # iterate over the layers
            for layer, params in layers.items():
                # get ymin and ymax of the layer
                ymin, ymax = params['y_coord']

                # check if part of the slice in the current layer
                if ymax > y:
                    # check if top of the layer is higher than slice
                    if ymax > ym_top:
                        height[layer] = (ym_top, y)
                    else:
                        height[layer] = (ymax, y)

                        # set tracker to top of the current layer
                        y = ymax

        elif len(layers) == 1:
            # only one layer defined
            # get name of the layer
            name = list(layers.keys())[0]

            # add to dict
            height[name] = (ym_top, ym_low)

        else:
            # no layers have been defined, raise InputError
            raise InputError('No layers have been defined')

        return height

    def _get_intersect(self, x):
        """ Get the intersection between a slice and the circle

        method returns the y coordinate of the intersection point
        between the circle and the slice. Note that two y coordinates
        are possible, this method returns the one below the centre of
        the circle. This is, because the slip circle is always located
        below the centre of the circle

        Parameters
        ----------
        x : float
            x coordinate of the intersection

        Returns
        -------
        float
            y coordinate of the intersection, note that this is the
            lowest intersection point
        """
        # compute a, b and c for the quadratic formula
        a = 1
        b = -2*self.xy[1]
        c = self.xy[1]**2 - self.r**2 + (x-self.xy[0])**2

        # compute factor of the sqrt
        D = b**2 - 4*a*c

        # check value of D
        if D > 0:
            # two solutions possible, compute y coordinates
            y1 = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
            y2 = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)

            # return lowest value
            return np.min([y1, y2])

        elif D == 0 or D > -0.1**10:
            # one solutions
            return -b/(2*a)

        else:
            # negative value, imaginary solutions, raise error
            raise ValueError(
                'When making slices an imaginary solutions is encountered')

    def _angle(self, yc_inter, x1, y1, x2, y2):
        """ Compute the slip angle of a slice

        Method computes the slip angle, :math:`\\alpha_s`, of a slice
        with respect to the centre of the slice

        Parameters
        ----------
        yc_inter : float
            y coordinate of the intersection between the circle and a
            line vertical from the centre of the circle
        x1, y1 : float
            x and y coordinate of the lower left point of the slice
        x2, y2 : float
            x and y coordinate of the lower right point of the slice

        Returns
        -------
        float
            slip angle of the slice [rad]
        """
        # compute coordinate of the middle of the slice
        xm = x1 + 0.5*(x2-x1)
        ym = y1 + 0.5*(y2-y1)

        # check if xm is equal to xc
        if xm == self.xy[0]:
            # angle is zero
            return 0
        else:
            # compute chord from the line vertical from the middle of
            # the slice to the middle of the slice
            chord = np.sqrt((xm-self.xy[0])**2 + (ym-yc_inter)**2)

            # compute slip angle of the middle of the slice
            alpha_s = 2*np.arcsin(0.5*chord/self.r)

            # check if right or left from the centre of the circle
            if xm < self.xy[0]:
                # left of the centre, angle must be negative
                return -alpha_s
            elif xm > self.xy[0]:
                # right from the centre, positive angle
                return alpha_s

    def _make_figure(self, ax, show_slices):
        """ Method to generate a figure of the circle """
        # plot circle
        circle1 = plt.Circle(
            self.xy, self.r, color='darkgrey', ls='--', fill=False)
        ax.add_artist(circle1)

        # plot centre of the circle
        ax.plot(self.xy[0], self.xy[1], marker='.', color='darkgrey')

        # check if slices must be shown
        if self.slices is not None and show_slices:
            # iterate over the slices
            for slice, coord in self.slices.items():
                # plot the slice
                ax.plot(coord['x'], coord['y'], color='k')

        elif self.slices is not None and show_slices:
            # slices have not been generated but are requested
            user_warning(
                ('show_slices is True, but the slices have not yet been '
                 'generated'))

        return ax

    def make_slices(self, num_slices, point2, layers):
        """ Make slices

        Method to divide the given circle into slices

        Parameters
        ----------
        num_slices : int
            number of slices
        circle : dict
            dictionary with the centre and radius of the circle
        layers : dict
            dictionary with the layers
        """
        # check if number of slices is an integer
        if not isinstance(num_slices, int):
            raise InputError('Slices must be an int')

        # compute width of one slice
        width = point2[0]/num_slices

        # set dict to store the coordinates of the slices
        slices = {}

        # add first slice
        # compute intersection with circle (rightmost intersection)
        y_inter = self._get_intersect(x=width)

        # compute y coordinate of the intersection between a line from
        # the centre of the circle and the circle
        yc_inter = self._get_intersect(x=self.xy[0])

        # compute slip angle of the first slice
        alpha_s = self._angle(yc_inter, 0, 0, width, y_inter)

        # compute the height of the slice
        heights = self._compute_heights(
            0, y_inter, point2[1]*width/point2[0], 0, layers)

        # add to dict
        slices[0] = {
            'x': [0, width, width, 0],
            'y': [0, y_inter, point2[1]*width/point2[0], 0],
            'alpha_s': alpha_s,
            'h': heights}

        # set rightmost intersection as leftmost intersection for
        # the next slice
        y_inter_1 = y_inter

        # compute the coordinates of each slice
        for slice in range(1, num_slices):
            # compute coordinates of the slice
            x1 = width*slice
            x2 = x1 + width
            y1 = point2[1]*x1/point2[0]
            y2 = point2[1]*x2/point2[0]

            # get intersection with the circle
            y_inter_2 = self._get_intersect(x=x2)

            # compute slip angle of the slice
            alpha_s = self._angle(yc_inter, x1, y_inter_1, x2, y_inter_2)

            # compute the height of the slice
            heights = self._compute_heights(
                y_inter_1, y_inter_2, y2, y1, layers)

            # add coordinates and slip angle to dict
            slices[slice] = {
                'x': [x1, x2, x2, x1, x1],
                'y': [y_inter_1, y_inter_2, y2, y1, y_inter_1],
                'alpha_s': alpha_s,
                'h': heights}

            # set rightmost intersection as leftmost intersection for
            # the next slice
            y_inter_1 = y_inter_2

        # set slices as attribute
        self.slices = slices

    def plot(self, show_slices=False):
        """ Plot the circle

        Parameters
        ----------
        show_slices : bool, optional, default: False
            if the slices must be shown, only available after the
            computation as the slices are not generated before the
            computation. The default value is False, meaning that the
            slices will not be plotted.
        """
        # set additional spacing
        spacing = self.r/5
        # create figure
        fig, ax = plt.subplots()

        # add circle
        ax = self._make_figure(ax, show_slices)

        # set layout
        ax.set_ylim(self.xy[1]-self.r-spacing, self.xy[1]+self.r+spacing)
        ax.set_xlim(self.xy[0]-self.r-spacing, self.xy[0]+self.r+spacing)
        ax.set_xlabel('x coordinate [m]')
        ax.set_ylabel('y coordinate [m]')
        ax.set_title('Slip cricle')
        ax.grid()
        fig.gca().set_aspect('equal', adjustable='box')
        fig.tight_layout()
        plt.show()
