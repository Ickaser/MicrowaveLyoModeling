# NOTE: Although this generates a result, it is not very rigorous.
# I have to assume things like Rp for a trehalose-mAb solution,
# as well as chamber pressure when only Pirani is measured.
# What it does do is demonstrate that figure 4 from Park et al. 2021 is way off.

using Roots
const LP = LyoPronto
using Latexify, UnitfulLatexify
using NonlinearSolve

plot_defaults_lprf()

# ----------------------------------------
# Load data into memory
MWdat = map(row->(t=row.t*u"hr", P=row.P*u"W"), CSV.read(datadir("exp_raw", "gitter2019_MW1b.csv"), Table))
T1dat = map(row->(t=row.t*u"hr", T=row.T*u"°C"), CSV.read(datadir("exp_raw", "gitter2019_T1b.csv"), Table))
pdat = map(row->(t=row.t*u"hr", pch=row.pch*u"μbar"|>u"mTorr"), CSV.read(datadir("exp_raw", "gitter2019_pch1b.csv"), Table))

fitdat = PrimaryDryFit(T1dat.t, T1dat.T .|> u"K", maximum(T1dat.T)|>u"K")
plot(fitdat, nmarks=40)

pch = RampedVariable(minimum(pdat.pch))
plot(pch)
@df pdat plot!(:t, :pch, ylim=(0, NaN))




# Get M1 params, for convenience
M1 = load(datadir("exp_pro", "M1_KvRpRf.jld2"))

begin
# --------------------------
# - Set up parameters
# Vial parameters
vialsize = "10R"
rad_i, rad_o = get_vial_radii(vialsize)
A_p = π*rad_i^2  # cross-sectional area inside the vial
A_v = π*rad_o^2 # vial bottom area
m_v = get_vial_mass(vialsize)
# Formulation and fill
c_solid = 0.15u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 1.4u"cm^2*hr*Torr/g" # Total guess
A1 = 16.0u"cm*hr*Torr/g" # Total guess
A2 = 1.0u"1/cm" # Total guess
Rp = RpFormFit(R0, A1, A2)
Vfill = 2.3u"mL"
# Heat transfer: ambient chamber wall with some radiation
K_shf_f = ConstPhysProp(1u"W/m^2/K") # Total guess
Tsh = RampedVariable(243.15u"K") # Assumed ambient chamber wall # Total guess
# Geometry
h_f0 = Vfill/A_p
m_f0 = Vfill * ρ_solution
# RF fit parameters, taken from M1
@unpack Bf, Bvw, Kvwf = M1
# Controllable inputs
f_RF = 2.45u"GHz"
P_per_vial = LinearInterpolation(MWdat.P /30, MWdat.t) # Have to guess how many vials. Tried 30

params_base = ParamObjRF((
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf_f, A_v, A_p),
    (pch, Tsh, P_per_vial),
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl),
    (f_RF, LyoPronto.eppf, LyoPronto.epp_gl),
    (Kvwf, Bf, Bvw),
))
end


u0 = ustrip.([u"g", u"K", u"K"], [m_f0, fitdat.Tfs[1][1], fitdat.Tfs[1][1]])
gensol = (x,tpf)->gen_sol_pd(x, tpf...; u0=u0)
p0 = [0.5, -2.3, -0.5]

trans_KBB = KBB_transform_basic(Kvwf, Bf, Bvw)
nls_g = (du,x,tpf)->LyoPronto.err_expT!(du, gensol(x, tpf), tpf[3])
err_nls = NonlinearFunction{true}(nls_g, resid_prototype=zeros(num_errs(fitdat)))
p0 = [0.0, -1.0, 0.1]
tsol = gensol(p0, (trans_KBB, params_base))
modrftplot(tsol, trimend=1)
plot!(fitdat)

# ub3_p = [10.0, 7.0, 7.0]
# lb3_p = -ub3_p

# @time objf_KBB(p0, (trans_KBB, params_base, fitdat))
# opt3 = solve(OptimizationProblem(objf_KBB, p0, (trans_KBB, params_base, fitdat)), optalg, maxiters=100, show_trace=true)
opt3 = solve(NonlinearLeastSquaresProblem(err_nls, p0, (trans_KBB, params_base, fitdat)), LevenbergMarquardt())
sol3 = gensol(opt3.u, (trans_KBB, params_base))
prm3 = sol3.prob.p


begin
pl = blankplot_hrC()
plot!(fitdat,nmarks=20, labsuffix=", experiment")
modrftplot!(sol3, trimend=1, sampmarks=false)
end

# -----------------
# Braatz group analytical model

const AM = LyoProntoNIIMBLRF.AnalyticalModel

# They assume a sublimation temperature of 256.15K
# This is equivalent to chamber pressure of 1300 μbar, since they do no mass transfer resistance
# I will instead use measured pressures to get temperature
Tm1 = find_zero( T->LyoPronto.calc_psub(T*u"K")-5u"μbar", 250)*u"K" # Lower bound on pch
Tm2 = find_zero( T->LyoPronto.calc_psub(T*u"K")-20u"μbar", 250)*u"K" # Upper bound on pch
Kv = 65u"W/m^2/K"  # from paper
Qppp_f = 242345u"W/m^3"  # from paper
kf = LyoPronto.k_ice
Tb0 = 236.85u"K" # from paper
Tbf = Tb0
r = 1u"K/minute"

pp1 = AM.Params(Kv, Qppp_f, kf, h_f0, r, Tb0, Tbf, Tm1)
pp2 = AM.Params(Kv, Qppp_f, kf, h_f0, r, Tb0, Tbf, Tm2)
t1, s, T1 = AM.calc_tsT(pp1);
t2, s, T2 = AM.calc_tsT(pp2);

begin
pl_all = blankplot_hrC()
@df T1dat exptfplot!(:t, :T, nmarks=40)
# @df thm_pd exptvwplot!(:t, :T3, nmarks=40)
modrftplot!(sol3, labsuffix=", LC-DIF",trimend=1)
plot!(t1, T1, c=:green, label="TLM, lower p")
plot!(t2, T2, c=:purple, label="TLM, upper p")
plot!(legend=:topleft)
# plot!([1, 2.7], [-37, -27], arrow=:arrow, c=:gray, linewidth=2, label="")
# pl_brtz = plot!(pl_all, u"hr", u"°C"; inset=bbox(0.38, 0.65, 0.3, 0.15), subplot=2)
# pl_brtz = pl_all[2]
# plot!(pl_brtz, t, T ,label="",  markersize=7, c=:green)
# plot!(pl_brtz; ylim=(-40.1, -39.4), xlabel="", ylabel="", )
end
savefig(plotsdir("gitter2019_compT.svg"))
savefig(plotsdir("gitter2019_compT.pdf"))
