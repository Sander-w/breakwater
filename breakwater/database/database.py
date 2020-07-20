import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from ..utils.exceptions import user_warning

try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    user_warning('The Basemap package is required to use BreakwaterDatabase') 

class BreakwaterDatabase:
    """ Breakwater Database

    Import the breakwater database developed by (Allsop et al., 2009)
    consist of completed breakwater projects with data ranging from
    design wave height to contractor. The constructed breakwaters are
    classified by breakwater type, the following types are currently
    included: Rubble Mound, Composite, Berm, Caisson and Revetments.

    The breakwater database is a separate module of :py:obj:`breakwater`
    and can be imported with the following command:

    .. code-block:: python

       from breakwater.database import BreakwaterDatabase

    .. note::
       To use this module the :py:obj:`Basemap` package is required,
       this dependency is additional to the dependencies mentioned in
       Section `2.1 <install.html#Dependencies>`__. See the following
       `link <https://matplotlib.org/basemap/users/installing.html>`__
       for an installation of :py:obj:`Basemap`

    Parameters
    ----------
    update : bool, optional, default: False
        if False the data is not updated and the included data is
        loaded, if True the database is loaded from the :py:attr:`source`

    Attributes
    ----------
    df : pd.DataFrame
        DataFrame of the breakwater database
    source : str
        url of the source
    """

    def __init__(self, update=False):
        """ See help(BreakwaterDatabase) for more info """
        # set source of the database as attribute
        self.source = 'http://kennisbank-waterbouw.tudelft.nl/breakwaters/printall.php'

        # check if data must be updated from the source
        if update:
            # updata data from the source
            dataframes = pd.read_html(self.source)

            # get the correct df from the dfs
            for dataframe in dataframes:
                # check df of the database by checking the number of columns
                if dataframe.shape[1] > 10:
                    # correct df
                    df = dataframe
                    break

            # check if headers are correct
            if df.columns.dtype != object:
                # set first row as header
                header = df.iloc[0]
                df = df[1:]
                df.columns = header

            # set id coluimn as index
            df.set_index('id', inplace=True)
            df.index = df.index.astype(np.int64)

            # fix column names of front and rear
            names = ['armour', '(unit)', 'size', 'slope']
            new_columns = []

            # iterate over the columns
            for i, column in enumerate(df.columns):
                # only columns after front and rear must be changed
                if (column == 'front' or column == 'rear'
                        and pd.isnull(df.columns[i+1:i+4]).any()):
                    # change names of columns
                    for name in names:
                        new_columns.append(f'{column} {name}')

                elif pd.isnull(column):
                    pass

                else:
                    new_columns.append(column)

            # update column names
            df.columns = new_columns

            # convert numeric values to numeric values
            df = df.apply(pd.to_numeric, errors='ignore')

            # format types
            types_fmt = {
                'Rubble mound': 'Rubble Mound',
                'rubble mound': 'Rubble Mound',
                'Composite breakwater': 'Composite',
                'composite': 'Composite',
                'berm breakwater': 'Berm',
                'Berm breakwater': 'Berm',
                'Berm Breakwater': 'Berm',
                'Caisson breakwater': 'Caisson',
                'caisson breakwater': 'Caisson',
                'caisson': 'Caisson'}
            df.type = df.type.replace(types_fmt)

            # fix incorrect classified types
            units = ['Tetrapod', 'Antifer', 'Cubes', 'Accropode',
                     'Accropode II', 'Core-Loc', 'Xbloc', 'XblocPlus',
                     'Dolos', 'COB', 'Stabit', 'Cubipod']
            for unit in units:
                df.type = df.type.replace(unit, 'Rubble Mound')

        else:
            # load data from csv file
            file = os.path.dirname(os.path.abspath(__file__))
            df = pd.read_csv(f'{file}\\database.csv', index_col=0)

        # set df as attribute
        self.df = df

        # set markers and colors as private attributes
        self._markers = {
            'Rubble Mound': 'o', 'Caisson': 's', 'Berm': 'D',
            'Revetment': '>', 'Composite': '*', 'Unclassified': 'p'}
        self._colors = {
            'Rubble Mound': '#1f77b4', 'Caisson': '#ff7f0e', 'Berm': '#2ca02c',
            'Revetment': '#d62728', 'Composite': '#9467bd',
            'Unclassified': '#bcbd22'}

    @staticmethod
    def _del_zeros(df):
        """ Delete rows containing a zero """
        # make list to store ids to delete
        to_delete = []

        # iterate over the df
        for id, row in df.iterrows():
            # iterate over the row
            for param, value in row.iteritems():
                # check if value is zero
                if value == 0:
                    # row must be deleted
                    to_delete.append(id)
                    break

        # delete rows and return df
        return df.drop(labels=to_delete, axis=0)

    @staticmethod
    def _show_unclassified(show, df):
        """ Helper method to filter the df for unclassified bw types """
        # check if unclassified must be plotted
        if show:
            # replace nan with unclassified so label is correct in plot
            df.type = df.type.replace(np.nan, 'Unclassified')
        else:
            # remove unclassified breakwaters
            df.dropna(subset=['type'], inplace=True)

        return df

    @staticmethod
    def _validate_excludes(input):
        """ Verify the input for the correct type """
        # check input of excludes
        if input is not None and not isinstance(input, list):
            # raise TypeError for incorrect type of exclude
            raise TypeError(
                ('Types to exclude must be given in a list, not as '
                f'{type(input).__name__}'))

    @staticmethod
    def _validate_all_excludes(input):
        """ Validate if specified types have not been excluded """
        if input:
            # not all given types were in the database,
            # and were thus not excluded, therefore show warning
            not_excluded = ', '.join(input)
            user_warning(
                ('The following types in exclude are not in the database: '
                f'{not_excluded}'))

    @property
    def unclassified(self):
        """ Get the number of unclassified breakwaters """
        # get the unclassified breakwaters
        return self.df.type.isna().sum()

    def report(self, save=True, save_path='data_report.xlsx', decimals=3):
        """ Make a report of the data in the database

        Method to generate a data report of the database. For each
        breakwater type the total and missing number of datapoints is
        determined. Furthermore, for numerical values the mean and
        standard deviation is also computed.

        .. note:
           Note that the mean and standard deviation is not computed
           for the following parameters: start, finish, front size,
           rear size, latitude and longitude

        Parameters
        ----------
        save : bool, optional, default: True
            If True an Excel version of the data report is generated,
            use :py:obj:`save_path` to specify a save path
        save_path : str
            File path to save the Excel file to
        decimals : int, optional, default: 3
            Number of decimal places to round to (default: 0). If
            decimals is negative, it specifies the number of
            positions to the left of the decimal point.

        Returns
        -------
        pd.DataFrame
            if :py:obj:`save` is False a DataFrame of the report is
            returned
        """
        # create df to store data report in
        data_report = pd.DataFrame()

        # drop unclassified breakwaters
        df = self._show_unclassified(False, self.df)

        # define columns of which the mean and standard deviation does have
        # to be computed
        no_computations_required = [
            'start', 'finish', 'front size', 'rear size', 'Lat', 'Lon']

        # make list for headers and subheaders
        headers, subheaders = ['Type'], ['']

        # set bool to track if subheaders have been added for wall column
        add_subheaders_wall = True

        # iterate over the bw types
        for i, bw_type in enumerate(df.type.unique()):
            # get only the bws of the current type
            bws = df[df.type == bw_type]

            # create dict with bw_type to store data
            type_data = {'type': [bw_type]}

            # iterate over the columns with dtype (for computing mean, etc)
            for column, dtype in bws.dtypes.iteritems():
                # add name of the column as a header if first iteration
                if i == 0 and column != 'type':
                    headers.append(column)

                # check if column has the bw_type
                if column == 'type':
                    # pass as this has already been added
                    continue

                # check type of the column
                elif dtype == np.object:
                    # column contains strings or is a coordinate
                    # count total number of datapoints
                    total = len(bws[column].values)

                    # check if column contains the slope
                    if 'slope' in column:
                        # check bw type
                        if bw_type == 'Caisson' or bw_type == 'Composite':
                            # have slope 1:0 but this is not needed since
                            # these are vertical structures
                            type_data[column] = ['-']

                            # add subheader if first iteration
                            if i == 0:
                                subheaders.append('')
                            continue

                        else:
                            # slope 1:0 is specified when the slope unknown
                            # so replace these values with nan
                            column_data = bws[column].replace('1:0', np.nan)

                            missing = column_data.isna().sum()

                    else:
                        # count missing datapoints
                        missing = bws[column].isna().sum()

                    # add to column
                    type_data[column] = [f'{total-missing}/{total}']

                    # add subheader if first iteration
                    if i == 0:
                        subheaders.append('')

                elif column in no_computations_required:
                    # replace zeros for nan
                    column_data = bws[column].replace(0, np.nan)

                    # compute total and missing datapoints
                    total = len(bws[column].values)
                    missing = column_data.isna().sum()

                    # add to column
                    type_data[column] = [f'{total-missing}/{total}']

                    # add subheader if first iteration
                    if i == 0:
                        subheaders.append('')

                else:
                    # column is an int or float
                    # check if bw type is not caisson and column is wall
                    if column == 'wall' and bw_type != 'Caisson':
                        # wall is not a parameter of other structures
                        type_data[f'{column} no'] = ['-']
                        type_data[f'{column} comp'] = ['-']

                        # add subheader if first iteration
                        if add_subheaders_wall:
                            headers.append(column)
                            subheaders.append('')
                            subheaders.append('')
                            add_subheaders_wall = False
                        continue

                    # replace zeros for nan
                    column_data = bws[column].replace(0, np.nan)

                    # compute total and missing datapoints
                    total = len(bws[column].values)
                    missing = column_data.isna().sum()

                    # compute average and standard deviation
                    mean = column_data.mean()
                    std = column_data.std()

                    # add data to dict
                    type_data[f'{column} no'] = [f'{total-missing}/{total}']
                    type_data[f'{column} comp'] = [f'{mean}±{std}']

                    # append extra column to headers and set subheaders
                    if i == 0:
                        headers.append(column)
                        subheaders.append('datapoints')
                        subheaders.append('μ ± σ')

                        # check if column is wall for setting the bool
                        if column == 'wall':
                            # set to False because has been added here
                            add_subheaders_wall = False

            # add dict to df
            temp_df = pd.DataFrame(data=type_data)
            data_report = data_report.append(
                temp_df, ignore_index=True, sort=False)

        # update columns of the df with double row of headers
        data_report.columns = [headers, subheaders]

        # check if df must be saved
        if save:
            # save df
            data_report.to_excel(save_path)
        else:
            # return df
            return data_report

    def correlation(self, param1, param2, bw_type=None, method='pearson'):
        """ Compute the correlation between two parameters

        Parameters
        ----------
        param1 : str
            name of parameter 1
        param2 : str
            name of parameter 2
        bw_type : str, optional, default: None
            if specified only the values of the given bw_type are
            considered
        method : {pearson, spearman}, optional, default: pearson
            method of correlation

        Returns
        -------
        tuple
            correlation between the two parameters and the p-value
        """
        # check if bw_type is given
        if bw_type is not None:
            # only for bw_type
            pass

        else:
            # for all bw types
            # get the two columns and remove zeros
            filtered_df = self._del_zeros(self.df[[param1, param2]])

        # method of correlation
        if method == 'pearson':
            # pearson method
            corr = stats.pearsonr(
                filtered_df[param1].values, filtered_df[param2].values)

        elif method == 'spearman':
            # spearman method
            corr = stats.spearmanr(
                filtered_df[param1].values, filtered_df[param2].values)

        else:
            raise NotImplementedError(
                f'{method} is not supported, must be pearson or spearman')

        return corr

    def cross_section(self, id, B=None, Rc=None, h=None, slope=None):
        """ Plot a cross section of the breakwater

        Method to plot a cross section of a breakwater in the database.
        The breakwater is selected by the id of the breakwater. In case
        data is missing to plot the breakwater it is possible to specify
        these as arguments.

        .. warning::
           plot function currently only supports Rubble Mound and
           caisson breakwaters

        Parameters
        ----------
        id : int
            id of the breakwater to plot
        B : float, optional, default: None
            specify custom crest width
        Rc : float, optional, default: None
            specify custom crest height
        h : float, optional, default: None
            specify a custom water level
        slope : tuple, optional, default: None
            specify a custom slope, must be specified as a tuple (V, H)

        Raises
        ------
        ValueError
            If data required to plot a cross section is missing
        """
        # set custom slope if not specified for protecting _validate
        if slope is None:
            slope = (0,0)

        # get the data of the bw
        bw = self.df[self.df.index == id]
        type = bw.type.values[0]

        # get hydraulic data
        depth = _validate_plot_vals('depth', bw['depth'].values[0], h)

        # get geometric data
        freeboard = _validate_plot_vals('Rc', bw['heigth'].values[0], Rc)
        width = _validate_plot_vals('B', bw['width'].values[0], B)

        # set infobox with general info
        cost = bw['cost(M$)'].values[0]
        infobox = '\n'.join((
            f'constructed between {bw.start.values[0]} and {bw.finish.values[0]}',
            f'cost = {cost} M$',
            f'owner = {bw.owner.values[0]}',
            f'contractor = {bw.contractor.values[0]}',
            f'consultant = {bw.consultant.values[0]}',
            f'Hs = {bw.Hs.values[0]} m',
            f'Tz = {bw.Tz.values[0]} s',
            f'Tp = {bw.Tp.values[0]} s',
            f'Rc = {freeboard} m',
            f'h = {depth} m',))

        # check type
        if type == 'Rubble Mound':
            # get the slope of the bw
            slope_databasse = bw['front slope'].values[0].split(':')
            V = _validate_plot_vals(
                'slope', float(slope_databasse[0]), slope[0])
            H = _validate_plot_vals(
                'slope', float(slope_databasse[1]), slope[1])

            # compute coordinates
            x1 = H*(depth+freeboard)/V
            xwlev = H*freeboard/V

            x = [-x1-0.5*width, -0.5*width, 0.5*width, x1+0.5*width]
            y = [0, depth+freeboard, depth+freeboard, 0]

            # set additional width for xmin and xmax
            increase_x = 0

            # add extra info of rubble mound
            armour = bw['front armour'].values[0]
            size = bw['front size'].values[0]
            unit = bw['front (unit)'].values[0]

            infobox = '\n'.join((
                infobox,
                f'armour = {armour} of {size} {unit}',
                f'slope = {V}:{H}'))

        elif type == 'Caisson':
            # compute coordinates
            xwlev = 0

            x = [-0.5*width, -0.5*width, 0.5*width, 0.5*width]
            y = [0, depth+freeboard, depth+freeboard, 0]

            # set additional width for xmin and xmax
            increase_x = 2*width

        else:
            raise NotImplementedError(
                f'{type} is currently not supported for plotting')

        # create figure
        fig, ax = plt.subplots(figsize=(10,5))

        # plot bw
        ax.plot(x,y, color='k', lw=1)

        # get the xmin and xmax
        xmin = ax.get_xlim()[0]*1.2 - increase_x
        xmax = ax.get_xlim()[1]*1.2 + increase_x

        # plot bottom and wlev (left and right)
        ax.axhline(y=0, color='peru', lw=1, zorder=5)
        ax.hlines(
            y=depth, xmin=xmin, xmax=-xwlev-0.5*width, color='dodgerblue', lw=1)
        ax.hlines(
            y=depth, xmin=xwlev+0.5*width, xmax=xmax, color='dodgerblue', lw=1)

        # place the info box with all info
        props = dict(boxstyle='round', facecolor='whitesmoke', alpha=0.5)
        ax.text(
            1.02, 0.99, infobox, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props)


        # format the figure
        ax.set_xlim(xmin, xmax)
        plt.title(
            (f'Cross section of the {type.lower()} breakwater at '
             f'{bw.harbour.values[0]}, {bw.country.values[0]}'))

        ax.set_aspect('equal', adjustable='box')
        ax.grid()
        fig.tight_layout()
        plt.show()

    def map(
            self, area=[], resolution='c', show_unclassified=False,
            exclude=None):
        """ Plot the breakwaters on a world map

        Method to plot all breakwaters with coordinates on a map of
        the world, or part of the world if an area is specified. Method
        uses :py:obj:`Basemap` to generate the map.

        Parameters
        ----------
        area : list, optional, default: []
            specify the coordinates of the area to plot. Use following
            format [llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat]
        resolution : str, optional, default: c
            resolution of the map to use. Can be c (crude), l (low),
            i (intermediate), h (high), f (full).
        show_unclassified : bool, optional, default: False
            True is unclassified breakwaters with a coordinate must be
            plotted, False is unclassified breakwaters must not be
            plotted.
        exclude : list, optional, default: None
            list of breakwater types to exclude from the plot
        """
        # validate exclude input for type
        self._validate_excludes(exclude)

        # create the figure
        plt.figure(figsize=(16,12))

        # check if an area is specified
        if not area:
            # plot full map
            m = Basemap(
                lat_0=0, lon_0=0, projection='robin', resolution=resolution)
            size = 7

        else:
            # plot the specified area
            m = Basemap(
                llcrnrlon=area[0], llcrnrlat=area[1], urcrnrlon=area[2],
                urcrnrlat=area[3], resolution=resolution)
            size = 15

        # edit lay-out of the map
        m.drawcountries(color='#d7d7d7')
        m.drawmapboundary(fill_color='#D0CFD4', linewidth=0)
        m.fillcontinents(color='#EFEFEF', lake_color='#D0CFD4')

        # remove bw's without a Lon and\or Lat
        filtered_df = self._del_zeros(self.df[['type', 'Lon', 'Lat']])

        # filter df for showing unclassified bw types
        filtered_df = self._show_unclassified(show_unclassified, filtered_df)

        # iterate over the bw types to plot them
        for i, bw_type in enumerate(filtered_df.type.unique()):
            # check if type is in exclude
            if exclude is not None and bw_type in exclude:
                # pass and delete type from exclude
                exclude.remove(bw_type)

            else:
                # get only the bws of the current type
                bws = filtered_df[filtered_df.type == bw_type]

                # plot the bws
                xpt, ypt = m(bws.Lon.values, bws.Lat.values)
                m.scatter(
                    xpt, ypt, s=size, alpha=1, label=bw_type,
                    c=self._colors.get(bw_type, '#17becf'),
                    marker=self._markers.get(bw_type, 'H'), zorder=2)

        # check if all types in exclude are in the database
        self._validate_all_excludes(exclude)

        # add legend, set tight_layout and show plot
        plt.legend(loc=1)
        plt.tight_layout()
        plt.show()

    def scatter(
            self, param1, param2, show_unclassified=False, exclude=None,
            min_data=5, xmax=None, ymax=None, bins_param1=10, bins_param2=10):
        """ Make a scatter plot of two parameters

        Method to generate a scatter plot with histograms for two
        parameters in the database.

        Parameters
        ----------
        param1 : str
            name of parameter 1
        param2 : str
            name of parameter 2
        show_unclassified : bool, optional, default: False
            True is unclassified breakwaters with a coordinate must be
            plotted, False is unclassified breakwaters must not be
            plotted.
        exclude : list, optional, default: None
            list of breakwater types to exclude from the plot
        min_data : int, optional, default: 5
            minimum number of datapoints required for plotting, if the
            data for a bw type is less than this limit it will be skipped
        xmax : float, optional, default: None
            maximum x coordinate of the scatter plot, by default this
            limit is automatically determined
        ymax : float, optional, default: None
            maximum y coordinate of the scatter plot, by default this
            limit is automatically determined
        bins_param1 : str
            number of bins for param2
        bins_param2 : str
            number of bins for param2
        """
        # validate exclude input for type
        self._validate_excludes(exclude)

        # filter the df to remove zeros
        filtered_df = self._del_zeros(self.df[['type', param1, param2]])

        # filter df for showing unclassified bw types
        filtered_df = self._show_unclassified(show_unclassified, filtered_df)

        # create the figure
        fig = plt.figure(figsize=(12,9))
        gs = GridSpec(4,4)

        # create the scatter plot
        scatter_plot = fig.add_subplot(gs[1:4,0:3])

        # make the histograms of the two parameters
        top_hist = fig.add_subplot(gs[0,0:3])
        right_hist = fig.add_subplot(gs[1:4,3])

        # create empty list for storing data of bw types
        hist_param1, hist_param2 = [], []

        # create empty list to store colors for histogram
        colors = []

        # iterate over the bw types
        for bw_type in filtered_df.type.unique():
            # check if type is in exclude
            if exclude is not None and bw_type in exclude:
                # pass and delete type from exclude
                exclude.remove(bw_type)

            else:
                # get only the bws of the current type
                bws = filtered_df[filtered_df.type == bw_type]

                # check if there is enough data
                if (len(bws[param1].values) <= min_data
                        or len(bws[param2].values) <= min_data):
                    # not enough data, show warning
                    user_warning(
                        f'{bw_type} was skipped because of a lack of data')

                else:
                    # enough data, plot on scatter
                    scatter_plot.scatter(
                        bws[param1].values, bws[param2].values, s=18,
                        c=self._colors.get(bw_type, '#17becf'), label=bw_type,
                        marker=self._markers.get(bw_type, 'H'), zorder=2)

                    # add used color to the list
                    colors.append(self._colors.get(bw_type, '#17becf'))

                    # add to hist lists
                    hist_param1.append(bws[param1].values)
                    hist_param2.append(bws[param2].values)

        # check if all types in exclude are in the database
        self._validate_all_excludes(exclude)

        # determine x and y limits, if not specified as arguments
        if xmax is None:
            # determine xmax
            xmax = np.round(scatter_plot.get_xlim()[1])

        if ymax is None:
            # determine ymax
            ymax = np.round(scatter_plot.get_ylim()[1])

        # generate the bins
        param1_bins = np.linspace(0, xmax, bins_param1)
        param2_bins = np.linspace(0, ymax, bins_param2)

        # add the histograms to the plot
        top_hist.hist(hist_param1, bins=param1_bins, color=colors)
        right_hist.hist(
            hist_param2, bins=bins_param2, orientation='horizontal',
            color=colors)

        # set x and y lims
        scatter_plot.set_xlim(0, xmax)
        scatter_plot.set_ylim(0, ymax)
        top_hist.set_xlim(0, xmax)
        right_hist.set_ylim(0, ymax)

        # add grid to all plots
        scatter_plot.grid()
        top_hist.grid()
        right_hist.grid()

        # set labels
        scatter_plot.set_xlabel(param1.capitalize())
        scatter_plot.set_ylabel(param2.capitalize())
        top_hist.set_ylabel('Frequency')
        right_hist.set_xlabel('Frequency')

        # remove ticks from histograms
        plt.setp(top_hist.get_xticklabels(), visible=False)
        plt.setp(right_hist.get_yticklabels(), visible=False)

        # add legend, set tight_layout and show plot
        scatter_plot.legend()
        plt.tight_layout()
        plt.show()

    def hist(
            self, param, show_unclassified=False, exclude=None, min_data=5,
            xmax=None, bins=10):
        """ Plot a histogram of a parameter

        Parameters
        ----------
        param : str
            name of the parameter
        show_unclassified : bool, optional, default: False
            True is unclassified breakwaters with a coordinate must be
            plotted, False is unclassified breakwaters must not be
            plotted.
        exclude : list, optional, default: None
            list of breakwater types to exclude from the plot
        min_data : int, optional, default: 5
            minimum number of datapoints required for plotting, if the
            data for a bw type is less than this limit it will be skipped
        xmax : float, optional, default: None
            maximum x coordinate of the histogram, by default this
            limit is automatically determined
        bins : int, optional, default: 10
            number of bins
        """
        # validate exclude input for type
        self._validate_excludes(exclude)

        # remove the zero values from the df
        df = self.df[['type', param]]
        df = df[df[param] != 0]

        # filter df for showing unclassified bw types
        filtered_df = self._show_unclassified(show_unclassified, df)

        # create lists to store labels, data and colors
        labels, data, colors = [], [], []

        # set variable to track maximum x
        xmax_computed = 0

        # iterate over the types in the df to get the data
        for bw_type in filtered_df.type.unique():
            # check if type is in exclude
            if exclude is not None and bw_type in exclude:
                # pass and delete type from exclude
                exclude.remove(bw_type)

            else:
                # get only the bws of the current type
                bws = filtered_df[filtered_df.type == bw_type]

                # check if there is enough data
                if len(bws[param].values) <= min_data:
                    # not enough data, show warning
                    user_warning(
                        f'{bw_type} was skipped because of a lack of data')

                else:
                    # enough data, plot on scatter
                    # add used color, data and label to the list
                    labels.append(bw_type)
                    data.append(bws[param].values)
                    colors.append(self._colors.get(bw_type, '#17becf'))

                    # get maximum value
                    if np.max(bws[param].values) >= xmax_computed:
                        xmax_computed = np.max(bws[param].values)

        # check if all types in exclude are in the database
        self._validate_all_excludes(exclude)

        # determine x limit, if not specified as arguments
        if xmax is None:
            # determine xmax
            xmax = np.round(xmax_computed)

        # generate the bins
        generated_bins = np.linspace(0, xmax, bins)

        # plot histogram
        plt.hist(data, bins=generated_bins, label=labels, color=colors)

        # format axis
        plt.xlim(0, xmax)
        plt.xlabel(f'{param}')
        plt.ylabel('Frequency')

        # set other lay-out
        plt.title(f'Frequency Histogram of {param}')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.show()


def _validate_plot_vals(param, val, given_val):
    """ Validate the values to plot """
    if val == 0:
        # check if val is given
        if given_val != 0:
            # return specified value
            return given_val
        else:
            # raise error
            raise ValueError(
                (f'No value for {param} in the database, use arguments to '
                  'specify custom value'))
    else:
        # return value
        return val
