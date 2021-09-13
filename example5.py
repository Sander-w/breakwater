import breakwater as bw
from Cubipod import Cubipod
"""
How to use the C02 Footprint functions for the Revetment at Energy Island
    Hydraulic conditions still need to be updated -> as well the gradings will probably be (Cupipod)
"""

battjes = bw.BattjesGroenendijk(Hm0= 14.4, h= 27.71, slope_foreshore=(1,100))
H2_per = battjes.get_Hp(0.02)

# define a limit state with hydraulic parameters, and the allowed damage
ULS = bw.LimitState(
    h= 27.71, Hm0= 14.4, H2_per=H2_per, Tp= 17.7, Tm= 13.3, T_m_min_1= 15,
    Sd= 5, Nod= 4, q= 10, label='ULS'
                    )

NEN = bw.RockGrading(rho= 2650) #standard gradings

cubipod = Cubipod()

# Design the Energy Island case with Cubipods.

configs = bw.Configurations(
    structure=['CRMR'], LimitState=ULS, rho_w=1025,
    slope_foreshore=(1,100), Grading=NEN, slope=(1, 8), B=(5, 8, 4),
    Dn50_core=(0.2, 0.4, 3), N=2100, ArmourUnit= cubipod, filter_rule = 'XblocPlus')

material_cost = {'type': 'Material', 'price': {'LMA_5/40': 5000, 'LMA_10/60': 6000, 'LMA_40/200': 70000, 'LMA_15/300': 8000, 'LMA_60/300': 9000,
              'HMA_300/1000': 10000, 'HMA_1000/3000': 11000, 'HMA_3000/6000': 12000, 'HMA_6000/10000': 13000,
              'HMA_10000/15000': 14000}, 'core_price': 400, 'unit_price': 500, 'concrete_price': 600, 'fill_price': 700,
            'transport_cost': 1000, 'Investment': None, 'length': None}

c02_cost = {'type': 'C02', 'price': {'LMA_5/40': 10, 'LMA_10/60': 20, 'LMA_40/200': 30, 'LMA_15/300': 40, 'LMA_60/300': 50,
              'HMA_300/1000': 60, 'HMA_1000/3000': 70, 'HMA_3000/6000': 80, 'HMA_6000/10000': 90,
              'HMA_10000/15000': 100}, 'core_price': 100, 'unit_price': 300, 'concrete_price': 50, 'fill_price': 30,
            'transport_cost': None, 'Investment': None, 'length': None}

cost_dicts = [material_cost, c02_cost]

for i in range(len(cost_dicts)):
    NEN.add_cost(type= cost_dicts[i]['type'], cost= cost_dicts[i]['price'])
    configs.add_cost(type = cost_dicts[i]['type'], core_price= cost_dicts[i]['core_price'], unit_price= cost_dicts[i]['unit_price'],
                     concrete_price= cost_dicts[i]['concrete_price'], fill_price= cost_dicts[i]['fill_price'],
                     transport_cost=cost_dicts[i]['transport_cost'], investment= cost_dicts[i]['Investment'],
                     length= cost_dicts[i]['length'])
    configs.cost_influence(cost_dicts[i]['type'])

save = False

if save:
    # check if, and which, warnings have been encountered
    configs.show_warnings()

    # export concepts to the design explorer
    configs.to_design_explorer(params=['c02_cost', 'material_cost', 'slope', 'class armour', 'B', 'Rc'])

    # save the configs to a .breakwaters file
    configs.to_breakwaters('example3')





