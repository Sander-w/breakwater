import os
import sys
import string
try:
    import xlsxwriter
except ImportError:
    pass

from ._kwarg_validator import _RM_vkwargs, _C_vkwargs
from .exceptions import InputError, NotSupportedError

def _parameters(workbook, structure):
    """ Make sheet with input for parameters """
    # set first row and col
    row = 0
    first_col = 0

    # parameters to skip as these are object and not parameters
    skip = ['LimitState', 'Grading', 'ArmourUnit', 'BermMaterial']

    # make sheet
    sheet1 = workbook.add_worksheet(name='Parameters')

    # get parameters for the given structure
    params = {}
    if 'RRM' in structure and 'CRM' in structure:
        params.update(_RM_vkwargs(type='both'))
    elif 'RRM' in structure:
        params.update(_RM_vkwargs(type='Rock'))
    elif 'CRM' in structure:
        params.update(_RM_vkwargs(type='ArmourUnit'))

    if 'RC' in structure and 'CC' in structure:
        params.update(_C_vkwargs(type='both'))
    elif 'RC' in structure:
        params.update(_C_vkwargs(type='Rock'))
    elif 'CC' in structure:
        params.update(_C_vkwargs(type='ArmourUnit'))

    # set format to use
    fmt_header = workbook.add_format(
        {'bold': True, 'locked': True, 'right': 1, 'bottom': 2,
         'valign': 'vcenter'})

    fmt_header_center = workbook.add_format(
        {'bold': True, 'align': 'center', 'locked': True, 'right': 1,
         'bottom': 2, 'valign': 'vcenter'})

    fmt_header_center_var = workbook.add_format(
        {'bold': True, 'align': 'center', 'locked': True, 'right': 1,
         'valign': 'vcenter'})

    column = workbook.add_format({'right': 1, 'locked': True})

    column_center = workbook.add_format(
        {'align': 'center', 'locked': False, 'right': 1})
    required_column  = workbook.add_format(
        {'bg_color': '#FFC7CE', 'font_color': '#9C0006', 'align': 'center',
         'locked': False, 'right': 1})

    non_varying = workbook.add_format(
        {'align': 'center', 'bg_color': '#A6A6A6', 'locked': True, 'right': 1})

    sheet1.set_column(first_col, first_col, 20)
    sheet1.set_column(first_col+1, first_col+2, 10)

    # Turn worksheet protection on.
    options = {
        'format_cells':          True,
        'select_locked_cells':   False,
    }
    sheet1.protect('', options)

    # set the header of the columns
    headers = ['Parameter', 'Value', 'Varying']
    fmt = [fmt_header, fmt_header_center, fmt_header_center_var]
    for col, header in enumerate(headers):
        letter = string.ascii_uppercase[first_col+col]
        if header == 'Varying':
            # horizontal merge
            letter2 = string.ascii_uppercase[first_col+col+2]
            merge_range = f'{letter}{row+1}:{letter2}{row+1}'
            index_varying = col
        else:
            # vertical merge
            merge_range = f'{letter}{row+1}:{letter}{row+2}'

        sheet1.merge_range(merge_range, header, fmt[col])

    # increment row for writing the next row
    row += 1

    # write subheaders of Varying
    sheet1.write_string(row, index_varying, 'Min', fmt_header_center)
    sheet1.write_string(row, index_varying+1, 'Max', fmt_header_center)
    sheet1.write_string(row, index_varying+2, 'Num', fmt_header_center)

    # add parameters to sheet
    row += 1
    for key, info in params.items():
        if key not in skip:
            # write parameter to the 1st column
            sheet1.write_string(row, first_col, key, column)

            # write value to 2nd column if not required
            if not info['Required']:
                # not a required parameter, so add default value
                default = info['Default']

                # check if default is a tuple since these are not supported
                if isinstance(default, tuple) or isinstance(default, list):
                    # convert to string
                    default = str(default)

                # write to cel
                sheet1.write(row, first_col+1, default, column_center)
            else:
                # required parameter, thus write empty
                sheet1.write_blank(row, first_col+1, None, column_center)

            # write varying parameters to third column
            form_constant = f'ISBLANK(B{row+1})'

            if not info['Constant']:
                # parameter is allowed to vary
                varying_format = column_center

                # check if requried to set the correct formula
                if info['Required']:
                    # add conditional format to required values allowed to vary
                    form_varying = f'SUMPRODUCT(--(C{row+1}:E{row+1}<>""))<3'
                    form = (f'IF(({form_varying})=TRUE, IF({form_constant}='
                             'TRUE, TRUE, FALSE), FALSE)')
                    sheet1.conditional_format(
                        row, first_col+1, row, first_col+1,
                        {'type': 'formula', 'criteria': form,
                            'format': required_column})

                # add formula for data validation
                form_validation = (f'IF((SUMPRODUCT(--(B{row+1}:E{row+1}<>"")'
                                    ')=4)=TRUE, "Parameter set as constant '
                                    'and varying parameter, please choose one'
                                    '.", "")')
                sheet1.write_formula(row, first_col+5, form_validation)

            else:
                # parameter is not allowed to vary
                varying_format = non_varying

                if info['Required']:
                    sheet1.conditional_format(
                        row, first_col+1, row, first_col+1,
                        {'type': 'formula', 'criteria': form_constant,
                        'format': required_column})

            sheet1.write_blank(row, index_varying, None, varying_format)
            sheet1.write_blank(row, index_varying+1, None, varying_format)
            sheet1.write_blank(row, index_varying+2, None, varying_format)

            # increment row
            row += 1

def _sheet_generator(workbook, name, headers, required=None, num_rows=25):
    """ Make sheet for with input

    Make input sheet for the LimitState, RockGrading or ArmourUnits
    """
    # set number of columns and start row
    num_cols = len(headers)
    first_col = 0
    row = 0

    # add sheet to workbook
    sheet1 = workbook.add_worksheet(name)

    # set cel formatting
    # headers
    fmt_header = workbook.add_format(
        {'bold': True, 'locked': True, 'right': 1, 'bottom': 2})

    fmt_centered_header = workbook.add_format(
        {'bold': True, 'align': 'center', 'locked': True, 'bottom': 2,
         'right': 1})

    fmt_centered_header_last = workbook.add_format(
        {'bold': True, 'align': 'center', 'locked': True, 'bottom': 2})

    # first column
    fmt_param = workbook.add_format(
        {'right': 1, 'locked': False, 'bottom': 1})

    fmt_param_locked = workbook.add_format(
        {'right': 1, 'locked': True, 'bottom': 1})

    fmt_param_last_row = workbook.add_format({'right': 1, 'locked': False})

    # other columns
    fmt_val = workbook.add_format(
        {'align': 'center', 'locked': False, 'bottom': 1, 'right': 1})

    fmt_val_last_col = workbook.add_format(
        {'align': 'center', 'locked': False, 'bottom': 1})

    fmt_val_last_row = workbook.add_format(
        {'align': 'center', 'locked': False, 'right': 1})

    fmt_val_last_row_col = workbook.add_format(
        {'align': 'center', 'locked': False})

    fmt_empty = workbook.add_format({'locked': False})

    # turn on worksheet protection
    options = {
        'format_cells':          True,
        'insert_rows':           True,
        'delete_rows':           True,
        'select_locked_cells':   True,
    }
    sheet1.protect('', options)

    # set width of the sheets and write headers
    for col, header in enumerate(headers):
        if col == 0:
            # first col has a different width and format
            sheet1.set_column(col+first_col, col+first_col, 20)
            sheet1.write_string(row, col+first_col, header, fmt_header)
        elif col == num_cols-1:
            # format last column
            sheet1.set_column(col+first_col, col+first_col, 11)
            sheet1.write_string(
                row, col+first_col, header, fmt_centered_header_last)
        else:
            sheet1.set_column(col+first_col, col+first_col, 11)
            sheet1.write_string(
                row, col+first_col, header, fmt_centered_header)

    row += 1

    # check if there are required parameters
    if required is not None:
        # set the required required parameter
        param_index = headers.index('Parameter')
        value_index = headers.index('Value')
        for param in required:
            sheet1.write_string(
                row, param_index+first_col, param, fmt_param_locked)
            sheet1.write_blank(
                row, value_index+first_col, None, fmt_val_last_col)

            # increment row
            row += 1

    # add empty rows for input
    for i in range(num_rows):
        # check if last row
        if i < num_rows-1:
            # not the last row thus regular format
            for col in range(num_cols):
                if col == 0:
                    # different format for first column
                    sheet1.write_blank(row, col+first_col, None, fmt_param)
                elif col == num_cols-1:
                    # >2 columns thus different format for last column
                    sheet1.write_blank(
                        row, col+first_col, None, fmt_val_last_col)
                else:
                    sheet1.write_blank(row, col+first_col, None, fmt_val)
        else:
            # different format for last row
            for col in range(num_cols):
                if col == 0:
                    # different format for first column
                    sheet1.write_blank(
                        row, col+first_col, None, fmt_param_last_row)
                elif col == num_cols-1:
                    # >2 columns thus different format for last column
                    sheet1.write_blank(
                        row, col+first_col, None, fmt_val_last_row_col)
                else:
                    sheet1.write_blank(
                        row, col+first_col, None, fmt_val_last_row)

        # increment row
        row += 1

def generate_excel(filepath, input='configurations', structure=None):
    """ Generate excel file for design input

    Parameters
    ----------
    filepath : str
        location to save the excel input file
    input : str, optional, default: configurations
        specify which type of input excel must be generated, possible
        arguments: configurations (default), parameters, LimitState,
        Grading, ArmourUnits
    structure : {RRM, CRM, RC, CC}, optional, default: None
        structure for which the input sheet must be generated. RRM
        for a rubble mound with rock as armour layer, CRM for a rubble
        mound with concrete armour units as armour layer, RC for a
        vertical (composite) breakwater with rock as armour layer for
        the foundation and CC for a vertical (composite) breakwater with
        concrete armour units as armour layer for the foundation.

    Raises
    ------
    TypeError
        if the extension of the filepath is not .xlsx
    """
    # check if xlsxwriter is imported
    if 'xlsxwriter' not in sys.modules:
        # module was not imported, raise error
        raise ModuleNotFoundError(
            'The module Xlsxwriter is a required dependency for this function')

    # check if input is in supported input
    supported_input = [
        'configurations', 'parameters', 'LimitState', 'Grading',
        'ArmourUnits']
    if input not in supported_input:
        possible = ', '.join(supported_input)
        raise NotSupportedError(
            f'Excel input for {input} cannot be generated, must be {possible}')

    # get the extension of the file
    extension = os.path.splitext(filepath)[1]

    # check if an extension is included
    if not extension:
        # if not add .xlsx to the filepath
        filepath = f'{filepath}.xlsx'
    else:
        # there is an extension
        if extension != '.xlsx':
            # invalid extension
            raise TypeError(f'Extension {extension} is invalid, must be .xlsx')
        else:
            pass

    # generate the workbook
    workbook = xlsxwriter.Workbook(filepath)

    # add sheets
    if 'parameters' in input or 'configurations' in input:
        # get specified structure input
        if structure is not None:
            # convert the input of structure to a list
            if isinstance(structure, list):
                # must be a list so no change
                structure = structure
            elif isinstance(structure, str):
                # convert single input to list
                structure = [structure]
        else:
            raise InputError(
                ('Argument: structure must be specified when input is '
                 'parameters or configurations'))

        _parameters(workbook, structure)

    if 'LimitState' in input or 'configurations' in input:
        _sheet_generator(
            workbook, name='LimitState', headers=['Parameter', 'Value'],
            required=['h', 'label'])

    if 'Grading' in input or 'configurations' in input:
        _sheet_generator(
            workbook, name='RockGrading',
            headers=['Rock Class', 'M50 Lower', 'M50 Upper', 'NLL', 'NUL'])

    if 'ArmourUnits' in input or 'configurations' in input:
        # check if structure is given
        if structure is not None:
            # determine name of the sheet
            if 'CRM' in structure and 'CC' in structure:
                names = ['ArmourUnit', 'BermMaterial']
            elif 'CRM' in structure:
                names = ['ArmourUnit']
            elif 'CC' in structure:
                names = ['BermMaterial']
            else:
                names = []
        else:
            names = ['ArmourUnit']

        for name in names:
            _sheet_generator(
                workbook, name=name, headers=['Volume', 'D', 'h', 'Vc'])

    # save and close workbook
    workbook.close()
