
plot_defaults_mlm()

# ----------------------------------------
# See section "Specific Surface Area Measurements" for these details

# Load data into memory

T1mw = map(row->(t=row.t*u"hr", T=row.T*u"°C"), CSV.read(datadir("exp_raw", "bhambhani2021_fig5", "T_a_A-Exterior.csv"), Table))
T2mw = map(row->(t=row.t*u"hr", T=row.T*u"°C"), CSV.read(datadir("exp_raw", "bhambhani2021_fig5", "T_a_A-Interior.csv"), Table))
T1c = map(row->(t=row.t*u"hr", T=row.T*u"°C"), CSV.read(datadir("exp_raw", "bhambhani2021_fig5", "T_b_TC1.csv"), Table))
T2c = map(row->(t=row.t*u"hr", T=row.T*u"°C"), CSV.read(datadir("exp_raw", "bhambhani2021_fig5", "T_b_TC5.csv"), Table))

Tsh = RampedVariable([-50.0, -50, -21]u"°C" .|>u"K", [0.5, 0.5]u"K/minute", [0.5u"hr"] ) # Conventional lyo temperature

blankplot_hrC()
plot!(Tsh, tmax=20u"hr", yunit=u"°C")
@df T1c plot!(:t, :T)
@df T2c plot!(:t, :T)

fitdat_c = PrimaryDryFit(T1c.t, (T1c.T[T1c.T .< -35u"°C"] .|> u"K",
    T2c.T[T2c.T .< -35u"°C"] .|> u"K"); t_end=18.0u"hr")
fitdat_mw = PrimaryDryFit(T1mw.t, (T1mw.T[T1mw.T .< -32u"°C"] .|> u"K",
    T2mw.T[T2mw.T .< -32u"°C"] .|> u"K"), 
    maximum(T2mw.T)|>u"K", (3.0u"hr", 4.5u"hr"))
plot(fitdat_c, nmarks=0, )
plot(fitdat_mw, nmarks=0, )

pch_c = RampedVariable(40u"mTorr")
pch_mw = RampedVariable(60u"mTorr")

# 152 vials in run for microwave
P_per_vial = ConstantInterpolation([400, 800, 1200, 800]u"W" ./152, 
    [0u"hr", 1u"hr", 3u"hr", 4u"hr"], dir=:left, extrapolation_right=ExtrapolationType.Constant)

# Get M1 params, for convenience
M1 = load(datadir("exp_pro", "M1_KvRpRf.jld2"))

begin
# --------------------------
# - Set up parameters
# Vial parameters
A_v = π*(17u"mm"/2)^2 # vial bottom area, listed as 17mm diameter in materials
A_p = π*(15u"mm"/2)^2  # cross-sectional area inside the vial: guessing 1mm thickness
m_v = 5.0u"g" # Approximately, as halfway between 2R and 4R
 
# Formulation and fill
c_solid = 0.05u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
Vfill = 0.7u"mL"
R0 = 1.4u"cm^2*hr*Torr/g" # Total guess
A1 = 16.0u"cm*hr*Torr/g" # Total guess
A2 = 1.0u"1/cm" # Total guess
Rp = RpFormFit(R0, A1, A2)

# Heat transfer in lyo: give a typical guess
K_shf_c = ConstPhysProp(12.0u"W/m^2/K") # Total guess
# Heat transfer in vacuum dryer: ambient chamber wall with some radiation
K_shf_mw = ConstPhysProp(1.0u"W/m^2/K") # Total guess
Tsh_mw = RampedVariable([203.15u"K", 283.15u"K"], 70u"K/minute") # Ambient chamber wall, starting at -70 for simulation convenience then going to 0

# Geometry
h_f0 = Vfill/A_p
m_f0 = Vfill * ρ_solution
# RF fit parameters, taken from M1 as guess
@unpack Bf, Bvw, Kvwf = M1
# Controllable inputs
f_RF = 2.45u"GHz"

po_conv = ParamObjPikal((
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf_c, A_v, A_p),
    (pch_c, Tsh),
))

po_mw = ParamObjRF((
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf_mw, A_v, A_p),
    (pch_mw, Tsh_mw, P_per_vial),
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl),
    (f_RF, LyoPronto.eppf, LyoPronto.epp_gl),
    (Kvwf, Bf, Bvw),
))
end

# ---------------
# First, get Rp and Kv for conventional

trans_KRp = KRp_transform_basic(K_shf_c(0), R0, A1, A2)
p0_c = [0.0, 1.1, 0.1, 0.1]


nls_c = NonlinearFunction{true, SciMLBase.FullSpecialize}(nls_pd!, resid_prototype=zeros(num_errs(fitdat_c)))
opt_c = solve(NonlinearLeastSquaresProblem(nls_c, p0_c, (trans_KRp, po_conv, fitdat_c)), LevenbergMarquardt())

inverse(trans_KRp, transform(trans_KRp, opt_c.u))
conv_prof = gen_sol_pd(opt_c.u, trans_KRp, po_conv)

begin
blankplot_hrC()
modconvtplot!(conv_prof, label="model")
# plot!(fitdat_c, nmarks=0, label="data")
exptfplot!(T1c.t, T1c.T, T2c.T, nmarks=0)
plot!(Tsh, tmax=20u"hr", label=L"T_{sh}", c=:black)
plot!(xlim=(0, 20), ylim=(-60, 40))
end

begin
plot(xunit="cm", yunit="cm^2*Torr*hr/g", ylabel=L"R_p", xlabel=L"h_d")
plot!(calc_hRp_T(po_conv, fitdat_c; i=1)..., label="TC1")
plot!(calc_hRp_T(po_conv, fitdat_c; i=2)..., label="TC6")
l = range(0u"cm", stop=h_f0, length=100)
plot!(l, conv_prof.prob.p.Rp.(l), label="fit to T")
end

Bi_z_conv = conv_prof.prob.p.Kshf.val * h_f0 / LyoPronto.k_ice |> NoUnits

# -----------------
# Microwave fitting

@reset po_mw.Rp = conv_prof.prob.p.Rp # Use Rp from conv. drying
trans_KBB = KBB_transform_bounded(Kvwf, Bf, Bvw)
trans_K = K_transform_basic(K_shf_mw(0))
trans_Rp = Rp_transform_basic(R0, A1, A2)
trans_KKBB = merge(trans_K, trans_KBB)
trans_KKBBRp = merge(trans_KKBB, trans_Rp)
nls_mw = NonlinearFunction{true, SciMLBase.FullSpecialize}(nls_pd!, resid_prototype=zeros(num_errs(fitdat_mw)))
p0_mw = [2.0, 1.0, -5.1, -1.5] # obj_expT 220.8
p0_mw = [-1.0, 0.0, -1.1, -2.0] # eventually: obj_expT 141
p0_mwRp = [2.0, 1.0, -0.1, -1.0, -1.0, 0.0, 4.0 ] 
tsol = gen_sol_pd(p0_mwRp, trans_KKBBRp, po_mw)
obj_expT(tsol, fitdat_mw, verbose=true)
modrftplot(tsol)
plot!(fitdat_mw)

opt_mw_sameRp = solve(NonlinearLeastSquaresProblem(nls_mw, p0_mw, (trans_KKBB, po_mw, fitdat_mw)), LevenbergMarquardt())
sol_mw_sameRp = gen_sol_pd(opt_mw_sameRp.u, trans_KKBB, po_mw)

opt_mw = solve(NonlinearLeastSquaresProblem(nls_mw, p0_mwRp, (trans_KKBBRp, po_mw, fitdat_mw)), LevenbergMarquardt())
sol_mw = gen_sol_pd(opt_mw.u, trans_KKBBRp, po_mw)

obj_expT(sol_mw, fitdat_mw, verbose=true)
obj_expT(sol_mw_sameRp, fitdat_mw, verbose=true)
begin
blankplot_hrC()
modrftplot!(sol_mw, label="mw Rp")
modrftplot!(sol_mw_sameRp, c=:purple, label="conv Rp")
plot!(fitdat_mw)
end
qplotrf(sol_mw)

transform(trans_KKBBRp, opt_mw.u)

Bi_z_rf = sol_mw.prob.p.Kshf.val * h_f0 / LyoPronto.k_ice |> NoUnits
Bi_r_rf = sol_mw.prob.p.Kvwf * 7.5u"mm" / LyoPronto.k_ice |> NoUnits
@show Bi_z_conv Bi_z_rf Bi_r_rf

save(datadir("exp_pro", "bhambhani2021_fit_params.jld2"), 
    Dict("conv"=>transform(trans_KRp, opt_c.u),
    # "mw"=>transform(trans_KKBB, opt_mw.u)))
    "mw"=>transform(trans_KKBBRp, opt_mw.u)))

Rp_conv = transform(trans_KRp, opt_c.u).Rp
Rp_mw = transform(trans_KKBBRp, opt_mw.u).Rp
l = range(0u"cm", stop=h_f0|>u"cm", length=101)

begin
plc = blankplot_hrC()
# plot!(fitdat_c, nmarks=0, label="data")
exptfplot!(T1c.t, T1c.T, T2c.T, nmarks=0)
modconvtplot!(conv_prof)
plot!(Tsh, tmax=20u"hr", label=L"T_\mathrm{sh}", c=:black)
# tendplot!(fitdat_c.t_end, ls=:dash, label="")
plot!(xlim=(0, 20), ylim=(-60, 45), title="Conventional, S1")
plmw = blankplot_hrC(ylim=(-60, 45))
# tendplot!(fitdat_mw.t_end, ls=:dash, label="")
exptfplot!(T1mw.t, T1mw.T, T2mw.T, nmarks=25)
# modrftplot!(sol_mw_sameRp, labsuffix=", conv \$R_p\$", c=:purple, alpha=0.5)
modrftplot!(sol_mw, )
modrftplot!(sol_mw_sameRp, labels=["\$T_\\mathrm{f}\$" "\$T_\\mathrm{vw}\$"] .*", conv. \$R_p\$", c=:purple, alpha=0.7, trimend=1)
# modrftplot!(sol_mw, labels=["fit \$R_p\$" ""], )
hline!([Tsh_mw.setpts[end]], label=L"T_\mathrm{sh}", c=:black)
plot!(title="Microwave, S2")
plot!(xlim=(0, 20), legend=:bottomright)
# plot!(fitdat_c, nmarks=0, label="data")
pl_q = qplotrf(sol_mw, ordering=[1,2,3])
# plot!(title="Microwave, S2")
plot!(ylabel="Vial Heating\n", xlim=(0, 5), ylim=(0, 0.25), legend=:topright)
# hline!([0], c=:gray, lw=1, label="")
pl_P = plot(xlim=(0,5), xunit=u"hr", yunit=u"W", xlabel="Time", ylabel="MW Per\nVial", unitformat=:square, )
plot!(P_per_vial, label="", fillrange=[0], lw=0, markeralpha=0, c=:red)
plot!(pl_q, xticks=(0:5, []), xlabel="", bottom_margin=-5Plots.px)
pl_Rp = plot(u"cm", u"cm^2*hr*Torr/g", ylabel="R_p", xlabel=L"h_d", yunitformat=(l,u)->"\$$l\$\n["*latexify(u)*"]", xunitformat=:square)
plot!(l, Rp_conv.(l), label="conv")
plot!(l, Rp_mw.(l), label="MW", legend=:right, bottom_margin=20Plots.px)
plot!(ylim=(0, 5), )
extrapls = plot(pl_q, pl_P, pl_Rp, layout=(3,1), left_margin=20Plots.px)


plot!(plc, left_margin=20*Plots.px)
plot!(plmw, ylabel="", left_margin=0*Plots.px, bottom_margin=20Plots.px)
plot(plc, plmw, extrapls, layout=(1,3), size=(1000,400))
end

savefig(plotsdir("bhambani2021_fit.svg"))
savefig(plotsdir("bhambani2021_fit.pdf"))


# ------------
# Version of plot for ISLFD East Coast 2025

begin
plc = blankplot_hrC()
# plot!(fitdat_c, nmarks=0, label="data")
exptfplot!(T1c.t, T1c.T, T2c.T, nmarks=0)
modconvtplot!(conv_prof)
plot!(Tsh, tmax=20u"hr", label=L"T_{sh}", c=:black)
# tendplot!(fitdat_c.t_end, ls=:dash, label="")
plot!(xlim=(0, 20), ylim=(-60, 45), title="Conventional")
plmw = blankplot_hrC(ylim=(-60, 45))
# tendplot!(fitdat_mw.t_end, ls=:dash, label="")
exptfplot!(T1mw.t, T1mw.T, T2mw.T, nmarks=25)
# modrftplot!(sol_mw_sameRp, labsuffix=", conv \$R_p\$", c=:purple, alpha=0.5)
modrftplot!(sol_mw, )
# modrftplot!(sol_mw_sameRp, labels=["\$R_p\$ of conv" ""], c=:purple, alpha=0.5)
# modrftplot!(sol_mw, labels=["fit \$R_p\$" ""], )
hline!([Tsh_mw.setpts[end]], label=L"T_{sh}", c=:black)
plot!(title="Microwave")
# plot!(fitdat_c, nmarks=0, label="data")
pl_q = qplotrf(sol_mw, ordering=[1,2,3])
# plot!(title="Microwave, S2")
plot!(ylabel="Vial Heating\n", xlim=(0, 5), ylim=(0, 0.25), legend=:topright)
# hline!([0], c=:gray, lw=1, label="")
pl_P = plot(xlim=(0,5), xunit=u"hr", yunit=u"W", xlabel="Time", ylabel="MW Per\nVial", unitformat=:square, )
plot!(P_per_vial, label="", fillrange=[0], lw=0, markeralpha=0, c=:red)
plot!(pl_q, xticks=(0:5, []), xlabel="", bottom_margin=-5Plots.px)
pl_Rp = plot(ylabel="R_p", xlabel=L"h_d", yunitformat=(l,u)->"\$$l\$\n["*latexify(u)*"]", xunitformat=:square)
plot!(l, Rp_conv.(l), label="conv")
plot!(l, Rp_mw.(l), label="MW", legend=:right, bottom_margin=20Plots.px)
plot!(ylim=(0, 5), )
extrapls = plot(pl_q, pl_P, pl_Rp, layout=(3,1), left_margin=20Plots.px)


plot!(plc, left_margin=20*Plots.px)
plot!(plmw, ylabel="", yticks=(-50:25:25, []), left_margin=-10*Plots.px, bottom_margin=20Plots.px)
plot(plc, plmw, extrapls, layout=(1,3), size=(1000,400))
end

savefig(plotsdir("bhambani2021_fit_slide.svg"))
savefig(plotsdir("bhambani2021_fit_slide.pdf"))
