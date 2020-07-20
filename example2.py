import breakwater as bw

# compute wave heights with Goda
H13_ULS, Hmax_ULS = bw.goda_wave_heights(
    h=15.1, d=12, Ho=5.3, T=9.4, slope_foreshore=(1,100))

H13_SLS, Hmax_SLS = bw.goda_wave_heights(
    h=12.1, d=9, Ho=3.3, T=7.9, slope_foreshore=(1,100))

# define ULS and SLS limit state with hydraulic parameters
ULS = bw.LimitState(
    h=15.1, H13=H13_ULS, Hmax=Hmax_ULS, T13=9.4, q=30, label='ULS')
SLS = bw.LimitState(
    h=12.1, H13=H13_SLS, Hmax=Hmax_SLS, T13=7.9, q=15, label='SLS')

# transform wave heights with relations from Rayleigh distribution
ULS.transform_periods(0.5)
SLS.transform_periods(0.5)

# define material
NEN = bw.RockGrading()

# design vertical breakwater with rock as foundation material
RC = bw.Caisson(
    Pc=0.2, rho_c=2400, rho_fill=1600, rho_w=1000, Bm=8, hb=2, layers=2,
    BermMaterial=NEN, LimitState=[ULS, SLS], slope_foreshore=(1,100), mu=0.5,
    beta=15)

# print the entire logger, to, for instance, see which formula is used
# to compute Rc
RC.print_logger(level='info')

# inspect the designed concept by plotting and printing all variants
RC.plot('all')
RC.print_variant('all')

# plot pressure distribution
RC.plot_pressure()
