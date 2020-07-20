import breakwater as bw

# compute wave heights not already known
battjes = bw.BattjesGroenendijk(Hm0=4.4, h=15, slope_foreshore=(1,100))
H2_per = battjes.get_Hp(0.02)

# define a limit state with hydraulic parameters, and the allowed damage
ULS = bw.LimitState(
    h=15, Hs=4.5, Hm0=4.4, H2_per=H2_per, Tp=9.4, Tm=8.8, T_m_min_1=9.7,
    Sd=5, Nod=2, q=20, label='ULS')

# define a material for the armour layer
NEN = bw.RockGrading(rho=2650)
xbloc = bw.Xbloc()

# design multiple configurations of RRM and CRM
configs = bw.Configurations(
    structure=['RRM', 'CRM'], LimitState=ULS, rho_w=1025,
    slope_foreshore=(1,100), Grading=NEN, slope=((1,3), (3,4), 4), B=(5, 8, 4),
    Dn50_core=(0.2, 0.4, 3), N=2100, ArmourUnit=xbloc)

# check if, and which, warnings have been encountered
configs.show_warnings()

# export concepts to the design explorer
configs.to_design_explorer(params=['slope', 'class armour', 'B', 'Rc'])

# save the configs to a .breakwaters file
configs.to_breakwaters('example3')
