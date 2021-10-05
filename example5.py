import breakwater as bw
from Cubipod import Cubipod

"""
How to use the C02 Footprint functions
"""
Hm0_315, Tp_315, Tm_315 = 14.2, 17.6, 13.5  # 1/10.000 yrs from 315 degrees North
Hm0_270, Tp_270, Tm_270 = 13.9, 16.4, 12.3  # 1/10.000 yrs from 270 degrees North
Hm0_225, Tp_225, Tm_225 = 12.7, 14.4, 11.1  # 1/10.000 yrs from 225 degrees North

battjes = bw.BattjesGroenendijk(Hm0=Hm0_225, h=27.71, slope_foreshore=(1, 100))
H2_per = battjes.get_Hp(0.02)

# define a limit state with hydraulic parameters, and the allowed damage
ULS = bw.LimitState(
    h=27.71,
    Hm0=Hm0_225,
    H2_per=H2_per,
    Tp=Tp_225,
    Tm=Tm_225,
    T_m_min_1=15,
    Sd=5,
    Nod=4,
    q=10,
    label="ULS",
)

NEN = bw.RockGrading(rho=2650)  # standard gradings

cubipod = Cubipod()

# Design the Energy Island case with Cubipods. Filter rule for single layer armoured Cubipod -> W/20

configs = bw.Configurations(
    structure=["CRMR"],
    LimitState=ULS,
    rho_w=1025,
    slope_foreshore=(1, 100),
    Grading=NEN,
    slope=(1, 3),
    B=(5, 8, 4),
    Dn50_core=(3, 4, 5),
    N=2100,
    ArmourUnit=cubipod,
    filter_rule="XblocPlus",
)

material_cost = {
    "type": "Material",
    "price": {
        "LMA_5/40": 5000,
        "LMA_10/60": 6000,
        "LMA_40/200": 70000,
        "LMA_15/300": 8000,
        "LMA_60/300": 9000,
        "HMA_300/1000": 10000,
        "HMA_1000/3000": 11000,
        "HMA_3000/6000": 12000,
        "HMA_6000/10000": 13000,
        "HMA_10000/15000": 14000,
    },
    "core_price": 400,
    "unit_price": 500,
    "concrete_price": 600,
    "fill_price": 700,
    "transport_cost": 1000,
    "Investment": None,
    "length": None,
}

cO2_cost = {
    "type": "CO2",
    "price": {
        "LMA_5/40": 10,
        "LMA_10/60": 20,
        "LMA_40/200": 30,
        "LMA_15/300": 40,
        "LMA_60/300": 50,
        "HMA_300/1000": 60,
        "HMA_1000/3000": 70,
        "HMA_3000/6000": 80,
        "HMA_6000/10000": 90,
        "HMA_10000/15000": 100,
    },
    "core_price": 100,
    "unit_price": 300,
    "concrete_price": 50,
    "fill_price": 30,
    "transport_cost": None,
    "Investment": None,
    "length": None,
}

cost_dicts = [material_cost, cO2_cost]

for i in range(len(cost_dicts)):
    NEN.add_cost(type=cost_dicts[i]["type"], cost=cost_dicts[i]["price"])
    configs.add_cost(
        type=cost_dicts[i]["type"],
        equipment= None,
        core_price=cost_dicts[i]["core_price"],
        unit_price=cost_dicts[i]["unit_price"],
        concrete_price=cost_dicts[i]["concrete_price"],
        fill_price=cost_dicts[i]["fill_price"],
        transport_cost=cost_dicts[i]["transport_cost"],
        investment=cost_dicts[i]["Investment"],
        length=cost_dicts[i]["length"],
    )
    configs.cost_influence(cost_dicts[i]["type"])

save = False

if save:
    # check if, and which, warnings have been encountered
    configs.show_warnings()

    # export concepts to the design explorer
    configs.to_design_explorer(
        params=["cO2_cost", "material_cost", "slope", "class armour", "B", "Rc"]
    )

    # save the configs to a .breakwaters file
    configs.to_breakwaters("example3")
