using DrWatson
@quickactivate :MicrowaveLyoModeling
# --------------- Set some plot defaults

plot_defaults_mlm()

# --------------- Read data
# Load processed data into memory
reloaded_SM1 = load(datadir("exp_pro", "SM1_processed.jld2"))
@unpack thm_conv_pd, lyo_conv_pd, fitdat_conv = reloaded_SM1

# ------ Model implementation ---------------------------------

begin
# Vial geometry
# rad = 1.0u"cm"
vialsize = "6R"
rad_i, rad_o = get_vial_radii(vialsize)
Ap = π*rad_i^2  # cross-sectional area inside the vial
Av = π*rad_o^2 # vial bottom area

# Formulation parameters
c_solid = 0.06u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 0.8u"cm^2*Torr*hr/g"
A1 = 14.0u"cm*Torr*hr/g"
A2 = 1.0u"1/cm"
Rpg = RpFormFit(R0, A1, A2)

# Cycle parameters
Vfill = 3u"mL" # ml
pch = RampedVariable(70u"mTorr")
T_shelf_0 = (273.15 -15 )u"K" # shelf temperature in K
T_shelf_final = (273.15 +10 )u"K"  # shelf temperature in K
ramp_rate = 0.5 *u"K/minute" # ramp rate deg/min
Tsh = RampedVariable([T_shelf_0, T_shelf_final], ramp_rate)

KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)

hf0 = Vfill / Ap

m_v = get_vial_mass(vialsize)
ρ_solution = 1u"g/mL"
m_f0 = Vfill * ρ_solution

f_RF = 8u"GHz"

# Guess values
K_vwf = 3.2e-3u"cal/s/K/cm^2"
B_f = 1.36e6u"Ω/m^2"
B_vw = 5.0e6u"Ω/m^2"

po_conv = ParamObjPikal((
    (Rpg, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh)
))
polc_conv = ParamObjRF((
    (Rpg, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh, RampedVariable(0u"W")), # dummy RF power
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl),
    (f_RF, ConstPhysProp(0.0), 0), # dummy frequency
    (12.0u"W/K/m^2", 1e6u"Ω/m^2", 1e6u"Ω/m^2") # dummy RF params
))
end

# -------------------

trans_Rp = Rp_transform_basic(R0, A1, A2)
u0 = ustrip.([u"g", u"K", u"K"], [m_f0, fitdat_conv.Tfs[1][1], fitdat_conv.Tfs[1][1]])
gensol = (x,tpf)->gen_sol_pd(x, tpf...; u0=u0)

p0_Rp = [1.1, 0.1, 0.1]
tsol1 = gensol(p0_Rp, (trans_Rp, polc_conv))
tsol2 = gen_sol_pd(p0_Rp, trans_Rp, po_conv)
begin
modconvtplot(tsol1, label="LC")
modconvtplot!(tsol2, label="Pikal", c=:green)
plot!(fitdat_conv, nmarks=10)
end

# opt_KRp = solve(OptimizationProblem(objf_KRp, p0_KRp, (po_conv, fitdat_conv)), optalg)
# conv_prof = gen_sol_KRp(opt_KRp.u, po_conv)

nls_SM1 = NonlinearFunction{true, SciMLBase.FullSpecialize}(nls_pd!, resid_prototype=zeros(num_errs(fitdat_conv)))
opt_Rp = solve(NonlinearLeastSquaresProblem(nls_SM1, p0_Rp, (trans_Rp, po_conv, fitdat_conv)), LevenbergMarquardt(), reltol=1e-8)
nls_SM1_lc = NonlinearFunction{true, SciMLBase.FullSpecialize}((du, x,tpf)->LyoPronto.err_expT!(du, gensol(x, tpf), tpf[3]), resid_prototype=zeros(num_errs(fitdat_conv)))
opt_Rp_lc = solve(NonlinearLeastSquaresProblem(nls_SM1_lc, p0_Rp, (trans_Rp, polc_conv, fitdat_conv)), LevenbergMarquardt())
# @time opt_Rp = solve(OptimizationProblem(objf_Rp, p0_Rp, (trans_Rp, porf_conv, fitdat_conv)), optalg)
# @show RpFormFit(transform(trans_Rp, opt_Rp.u).Rp...)
# @show RpFormFit(transform(trans_Rp, opt_Rp_lc.u).Rp...)
conv_prof_lc = gensol(opt_Rp_lc.u, (trans_Rp, polc_conv))
conv_prof = gen_sol_pd(opt_Rp.u, trans_Rp, po_conv)

begin
conv_labels = [L"T_{f1}"*", edge";;
               L"T_{f2}"*", exp.";;
               L"T_{f3}"*", exp.";;
]
blankplot_hrC()
@df thm_conv_pd exptfplot!(:t, :T4, :T1, :T2, linealpha=0.9, labels=conv_labels)
# plot!(fitdat_conv)
modconvtplot!(conv_prof)
modconvtplot!(conv_prof_lc, c=:green, label="LC")
plot!(Tsh, tmax=15u"hr", label=L"T_{sh}", c=:black)
plot!(xlim=(0, 14), ylim=(-35, 12), legend=:topright)
tendplot!(fitdat_conv.t_end)
end


# RF Data and Model ------------------------------------ 

reloaded_SM2 = load(datadir("exp_pro", "SM2_processed.jld2"))
@unpack thm_rf_pd, lyo_rf_pd, fitdat_rf = reloaded_SM2

# ---- Set up model for RF
# @info "Kv" Kshf_fit Kshf(pch(0))
Rp = transform(trans_Rp, opt_Rp.u).Rp


P_per_vial = RampedVariable(25u"W"/30*0.54)

po_rf = ParamObjRF((
    (Rp, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh, P_per_vial), 
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl),
    (f_RF, LyoPronto.eppf, LyoPronto.epp_gl),
    (K_vwf, B_f, B_vw),
))
# -----------------
@df lyo_rf_pd plot(:t, :Tsh_i)
plot!(Tsh, xlim=(0,2))

trans_KBB = KBB_transform_bounded(K_vwf, B_f, B_vw)

# gensolrf_SM = (x,tpf)->gen_sol_pd(x, tpf...; u0=u0_rf)
# obj_SM(x,p) = obj_expT(gensolrf_SM(x, p), p[end], tweight=1)
nls_SM2 = NonlinearFunction{true, SciMLBase.FullSpecialize}((du, x,tpf)->LyoPronto.err_expT!(du, gensol(x, tpf), tpf[3], tweight=1), resid_prototype=zeros(num_errs(fitdat_rf)))

p0_rf = [1, 5.0, -1.0]
tsol = gensol(p0_rf, (trans_KBB, po_rf))


opt2 = solve(NonlinearLeastSquaresProblem(nls_SM2, p0_rf, (trans_KBB, po_rf, fitdat_rf)), LevenbergMarquardt())

# opt1 = solve(OptimizationProblem(objf_SM, p0_rf, (po_rf, fitdat_rf)), optalg, maxiters=100, show_trace=true)
prof_RF = gensol(opt2.u, (trans_KBB, po_rf))

modrftplot(prof_RF)
plot!(fitdat_rf)

save_fitresults(opt2, "SM")

begin
convplot = blankplot_hrC()
@df thm_conv_pd exptfplot!(:t, :T4, :T2, :T1, nmarks=40, sampmarks=true, linealpha=0.2)

modrftplot!(conv_prof_lc, )
# plot!(edge_prof.t.*u"hr", edge_prof.(edge_prof.t, idxs=2).*u"K", label=L"T_{f}"*", model edge", c=6)
plot!(Tsh, label=L"T_\mathrm{sh}", c=:black)
plot!(xlim=(0, 12), ylim=(-35, 60), legend=:topleft)
tendplot!(fitdat_conv.t_end, ls=:dash, label="")

annotate!(7, 20, Plots.text("end of drying", 12, "Computer Modern"))
plot!([fitdat_conv.t_end-1u"hr", fitdat_conv.t_end-0.2u"hr"], [20, 20], arrow=:arrow, c=:black, linewidth=1, label="")
plot!(size=(400,400))
end

begin
rfplot = blankplot_hrC()
@df thm_rf_pd exptfplot!(:t, :T2, :T4, sampmarks=true, linealpha=0.2, nmarks=20)
@df thm_rf_pd exptvwplot!(:t, :T1, nmarks=20)
modrftplot!(prof_RF)
plot!(Tsh, c=:black, label=L"T_\mathrm{sh}")
plot!(size=(500,400), legend=:top)
tendplot!(fitdat_rf.t_end, ls=:dash, label="")
annotate!(9, 0, Plots.text("end of drying", 12, "Computer Modern"))
plot!([fitdat_rf.t_end+1.5u"hr", fitdat_rf.t_end+0.2u"hr"], [0, 0], arrow=:arrow, c=:black, linewidth=1, label="")
t_rf_off = thm_rf_pd.t[argmax(thm_rf_pd.T2)]
tendplot!(t_rf_off, ls=:dash, label="")
annotate!(8.0, -10, Plots.text("RF off", 12, "Computer Modern"))
plot!([t_rf_off+2.2u"hr", t_rf_off+0.2u"hr"], [-10, -10], arrow=:arrow, c=:black, linewidth=1, label="")
plot!(legend=:topright, ylim=(-35, 60), xlim=(0,12), size=(400,400))
end

# begin
# plot!(convplot, left_margin=20*Plots.px,legend=:none)
# plot!(rfplot, legend=(0.8, 0.7), ylabel="")
# plT = plot(convplot, rfplot, title = ["CONV" "RF"], layout=(1, 2),
#      size=(1000, 400), top_margin = 10Plots.px, bottom_margin=20Plots.px)
# end
# savefig(plotsdir("SM1-2_compT.svg"))
# savefig(plotsdir("SM1-2_compT.pdf"))

begin
plq = qplotrf(prof_RF, ordering=[3,2,1])
plot!(size=(400,300), ylim=(0, 0.4), widen=false)
end
# savefig(plotsdir("SM2_energy_budget.svg"))
# savefig(plotsdir("SM2_energy_budget.pdf"))

begin
plot!(convplot, left_margin=30*Plots.px,)
plot!(rfplot, ylabel=nothing, )
plot!(plq, left_margin=30Plots.px)
plot(convplot, rfplot, plq, title = ["SM1 (conv.)" "SM2 (RF-assisted)" "SM2"], layout=@layout([a b c{0.25w}]),
     size=(1100, 400), top_margin = 20Plots.px, bottom_margin=30Plots.px)
end
savefig(plotsdir("SM1-2_combine_T_q.svg"))
savefig(plotsdir("SM1-2_combine_T_q.pdf"))