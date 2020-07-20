# note running the entire script at once will result in an error,
# as the Excel file is first generated and must then be filled by
# the user. So comment out line 12 when running this script
# for the first time, and line 9 to prevent overwriting your input

import breakwater as bw

# generate the Excel file
bw.generate_excel('config input.xlsx', structure=['RRM', 'CRM'])

# load the excel file to make a design for multiple configurations
configs = bw.read_configurations(
    'config input.xlsx', structure=['RRM', 'CRM'], kd=16, name='Xbloc')
