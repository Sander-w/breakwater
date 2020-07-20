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

# design a conventional rubble mound breakwater with rock as armour layer
RRM = bw.RockRubbleMound(
    slope=(2,3), slope_foreshore=(1,100), rho_w=1025, B=5.5, N=2100,
    LimitState=ULS, Grading=NEN, Dn50_core=0.4)

# print the logger to see if any warnings were encountered during the design
RRM.print_logger(level='warnings')

# inspect the designed concept by plotting and printing all variants
RRM.plot('all')
RRM.print_variant('all')

# check validity of the used formulae
RRM.check_validity()

# define a concrete armour unit as armour layer
xbloc = bw.Xbloc()

# design a conventional rubble mound breakwater with xbloc as armour layer
CRM = bw.ConcreteRubbleMound(
    slope=(2,3), slope_foreshore=(1,100), B=5.5, rho_w=1025, LimitState=ULS,
    ArmourUnit=xbloc, Grading=NEN, Dn50_core=0.4)

# print the logger to see if any warnings were encountered during the design
CRM.print_logger(level='warnings')

# inspect the designed concept by plotting and printing all variants
CRM.plot('all')
CRM.print_variant('all')
