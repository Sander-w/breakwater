from ..utils.exceptions import user_warning

def scour_protection(L, slope=None):
    """ Compute the length of the scour protection

    Compute the estimated length of the scour hole with Sumer and
    Fredsoe (2000). In the table the estimated length of the scour hole
    from the experimental study is presented.

    +-------------------------+------------------+
    | slope of the breakwater | estimated length |
    +=========================+==================+
    | Vertical breakwater     | 1.0(L/4)         |
    +-------------------------+------------------+
    | 1:1.2                   | 0.6(L/4)         |
    +-------------------------+------------------+
    | 1:1.75                  | 0.4(L/4)         |
    +-------------------------+------------------+

    Parameters
    ----------
    L : float
        wave length at the toe of the structure, computed with the mean
        wave period [m]
    slope : float, optional, default: None
        Slope of the armour layer (V, H), for example a slope of 3V:4H
        is defined as (3, 4). By default the length of
        the scour protection for a vertical breakwater is computed.

    Returns
    -------
    float
        the required length of the scour protection
    """
    # check if a slope is given
    if slope is None:
        # no slope, thus vertical breakwater
        W = L/4

    else:
        # slope is given
        slope_div = slope[0]/slope[1]

        # determine scalar for estimating the length
        scalar = 0.4 + 0.2 * ((slope_div - 1/1.75)/(1/1.2 - 1/1.75))

        # compute length of the scour hole
        W = scalar*(L/4)

        if slope_div >= 1/1.75 and slope_div <= 1/1.2:
            # in range of experiments
            pass
        else:
            # slope is out of range, so print a warning
            user_warning(
                (f'slope {slope[0]}:{round(slope[1],2)} is out of range in '
                  'scour_protection(), must be between 1:1.75 and 1:1.2'))

    return W
