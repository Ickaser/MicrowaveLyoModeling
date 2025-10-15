plot_defaults_lprf()
# ----------------------------
# Read in data

data = load(datadir("exp_pro", "M3_processed.jld2"))
@unpack lyo_pd, thm_pd, thm_full, fitdat, P_per_vial, t_end = data
# Note: T4 and T3 frozen, T1 wall

# ---------------
# Set up simulation
begin
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)
# Mass transfer
R0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16u"cm*hr*Torr/g"
A2 = 0u"1/cm"
Rp = RpFormFit(R0, A1, A2)
# Other parameters
m_v = 7.6u"g"
c_solid = .05u"g/mL"
A_p = π* (1u"cm")^2
A_v = 3.8u"cm^2"
ρ_solution = 1u"g/mL"
# Geometry
Vfill = 3u"mL"
h_f0 = Vfill/A_p
m_f0 = Vfill * ρ_solution
# RF parameters
f_RF = 18u"GHz"

lcfit_M1 = load(datadir("exp_pro", "M1_KvRpRf.jld2"))
@unpack Bf, Bvw, Kvwf = lcfit_M1

pch = RampedVariable(100u"mTorr")
Tsh = RampedVariable(258.15u"K")
# P_per_vial loaded with data above

po = ParamObjRF((
    (Rp, h_f0, c_solid, ρ_solution),
    (Kshf, A_v, A_p),
    (pch, Tsh, P_per_vial), 
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl),
    (f_RF, LyoPronto.eppf, LyoPronto.epp_gl),
    (Kvwf, Bf, Bvw),
));
nothing
end

# -----------------

trans_KBB = KBB_transform_bounded(100u"W/m^2/K", Bf, Bvw)
u0 = ustrip.([u"g", u"K", u"K"], [m_f0, fitdat.Tfs[1][1], fitdat.Tfs[1][1]])
gensol = (x,tpf)->gen_sol_pd(x, tpf...; u0=u0)
p0 = [0.5, -2.3, -0.5]

@time sol = gensol(p0, (trans_KBB, po))

# Check the solution
modrftplot(sol)
plot!(fitdat)

# ---------------- Objective functions
# modrftplot(gensol(fill(1.0, 3), (po, fitdat)))
# obj_M3(x,tpf) = obj_expT(gensol(x, tpf), tpf[3], tweight=1)
# objf_M3 = OptimizationFunction(obj_M3, AutoForwardDiff(chunksize=3))
# @time opt = solve(OptimizationProblem(objf_M3, p0, (trans_KBB, po, fitdat)), optalg, )
nls_M3(du, x, tpf) = LyoPronto.err_expT!(du, gensol(x, tpf), tpf[3], tweight=1)
nlsf_M3 = NonlinearFunction{true}(nls_M3, resid_prototype=zeros(num_errs(fitdat)))
@time opt = solve(NonlinearLeastSquaresProblem(nlsf_M3, p0, (trans_KBB, po, fitdat)), LevenbergMarquardt())
prof_RF = gensol(opt.u, (trans_KBB, po, fitdat));

modrftplot(prof_RF)
# modrftplot!(prof_g, c=:red, labsuffix=" (guess)")
plot!(fitdat)

# Save fit results to a file
save_fitresults(opt, "M3")

begin
blankplot_hrC(legend=:bottomright)
modrftplot!(prof_RF, markeralpha=0)
@df thm_pd exptfplot!(:t, :T4, :T3)
@df thm_pd exptvwplot!(:t, :T1)
tendplot!(fitdat.t_end)
end



# -------------------- Nice plotting

begin
plT = blankplot_hrC()
plot!(Tsh, c=:black, label=L"T_\mathrm{sh}")
@df thm_pd exptfplot!(:t, :T4, :T3, nmarks=35)
@df thm_pd exptvwplot!(:t, :T1, nmarks=40, label=L"$T_\mathrm{vw1}$, exp.")
modrftplot!(prof_RF)
plot!(size=(480,400))
tendplot!(fitdat.t_end, label="", ls=:dash)
plot!(legend=:bottomleft, legend_columns=2, ylim=(-45, 30), xlim=(0,7.5), xlabel="", xformatter=x->"")

# annotate!(3.6, 20, Plots.text("end of drying", 12, "Computer Modern"))
# plot!([fitdat.t_end-2u"hr", fitdat.t_end-0.2u"hr"], [20, 20], arrow=:arrow, c=:black, linewidth=1, label="")
annotate!(5.8, -33, Plots.text("end of\n drying", 12, "Computer Modern"))
plot!([fitdat.t_end-1.2u"hr", fitdat.t_end-0.2u"hr"], [-40, -40], arrow=:arrow, c=:black, linewidth=1, label="")

plp = plot(u"hr", u"W", ylabel="RF Power", xlabel="Time", unitformat=:square)
plot!(P_per_vial.t, P_per_vial.u, c=:red, lw=0, fillrange=0, label="")
plot!(xlim=(0,7.5), ylim=(0, 3))
tendplot!(fitdat.t_end, label="", ls=:dash)
plot!(ylabel="RF Power \n[W/vial]", yunit=nothing, top_margin=-30Plots.px)
plTp = plot(plT, plp, layout=@layout([a; b{0.2h}]))
# end
# # savefig(plotsdir("M3_compT.svg"))
# # savefig(plotsdir("M3_compT.pdf"))

# begin
plq = qplotrf(prof_RF)
plot!(size=(400,300), xwiden=false, ylim=(0, 1.5), )
# savefig(plotsdir("M3_energy_budget.svg"))
# savefig(plotsdir("M3_energy_budget.pdf"))

plot(plTp, plq, layout=@layout([a  b{0.4w}]), size=(800, 400), bottom_margin=20Plots.px,left_margin=20Plots.px)
end
# plot(plT, plq, plp, layout=@layout([a  b{0.4w}; c{0.2h} d]), size=(1000, 400))
savefig(plotsdir("M3_T_q_combine.svg"))
savefig(plotsdir("M3_T_q_combine.pdf"))
