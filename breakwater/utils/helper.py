def round_even(value):
    """ Function to round values to the nearest even number

    Parameters
    ----------
    value : float
        value to round

    Returns
    -------
    int
        value rounded to the nearest even number
    """
    # round given value
    rounded_value = round(value)

    # check if number is even
    if (rounded_value % 2) == 0:
        # already an even number
        return rounded_value
    else:
        # value is an odd number
        # check if nearest even value is higher
        if (value - rounded_value) >= 0:
            # nearest even number is higher
            return int(rounded_value) + 1
        else:
            # nearest even number is lower
            return int(rounded_value) - 1
