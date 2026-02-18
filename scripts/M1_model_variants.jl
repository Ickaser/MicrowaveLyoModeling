using Roots

plot_defaults_mlm()

# ----------------------------------------
# Load processed data into memory
reloaded = load(datadir("exp_pro", "M1_processed.jld2"))
@unpack thm_pd, lyo_pd, fitdat = reloaded

begin
# --------------------------
# - Set up parameters
# Vial parameters
vialsize = "6R"
rad_i, rad_o = get_vial_radii(vialsize)
A_p = π*rad_i^2  # cross-sectional area inside the vial
A_v = π*rad_o^2 # vial bottom area
m_v = get_vial_mass(vialsize)
# Formulation and fill
c_solid = 0.05u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16.0u"cm*hr*Torr/g"
A2 = 0.0u"1/cm"
Rp = RpFormFit(R0, A1, A2)
Vfill = 5u"mL"
# Heat transfer
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
K_shf_f = RpFormFit(KC, KP, KD)
# Geometry
h_f0 = Vfill/A_p
m_f0 = Vfill * ρ_solution
# RF fit parameters (dummy values)
Bf = 2.0e7u"Ω/m^2"
Bvw = 0.9e7u"Ω/m^2"
Kvwf = 10.0u"W/K/m^2"
# Controllable inputs
f_RF = 8u"GHz"
pch = RampedVariable(100u"mTorr")
Tsh = RampedVariable([233.15u"K", 283.15u"K"], 0.5u"K/minute",)
P_per_vial = RampedVariable(10u"W"/17 * 0.54) # actual power / vial

params_base = ParamObjRF((
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf_f, A_v, A_p),
    (pch, Tsh, P_per_vial),
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl, A_v),
    (f_RF, LyoPronto.eppf, LyoPronto.epp_gl),
    (Kvwf, Bf, Bvw, missing),
))
end

prob = ODEProblem(params_base)
@time tsol = solve(prob, Rosenbrock23())

# -----------------
# Tune with model version 3 (DIF), no mass transfer params 


trans_KBB = KBB_transform_bounded(Kvwf, Bf, Bvw)
err_nls = NonlinearFunction{true}(nls_pd!, resid_prototype=zeros(num_errs(fitdat)))
p0 = [2.0, 3.8, 0.1]
tsol = gen_sol_pd(p0, trans_KBB, params_base)
modrftplot(tsol)
plot!(fitdat)

# ub3_p = [10.0, 7.0, 7.0]
# lb3_p = -ub3_p

# @time objf_KBB(p0, (trans_KBB, params_base, fitdat))
# opt3 = solve(OptimizationProblem(objf_KBB, p0, (trans_KBB, params_base, fitdat)), optalg, maxiters=100, show_trace=true)
opt3 = solve(NonlinearLeastSquaresProblem(err_nls, p0, (trans_KBB, params_base, fitdat)), LevenbergMarquardt())
sol3 = gen_sol_pd(opt3.u, trans_KBB, params_base)
prm3 = sol3.prob.p

# save_fitresults(opt3, "M1")

begin
blankplot_hrC()
modrftplot!(sol3)
plot!(fitdat)
end

#-------------------

# Tune with same model, but also fitting  a1-a2

aa_trans = as((;Rp=as((R0=Constant(Rp.R0), A1=TVScale(Rp.A1) ∘ TVExp(), A2=TVScale(0.01u"cm^-1") ∘ TVExp()))))
trans_KBBaa = merge(trans_KBB, aa_trans)

p0_3b = vcat(opt3.u, [0.0, -1.0])
# ub3b = [10.0, 7.0, 7.0, 3.0, 3.0]
# lb3b = -ub3b
opt3b = solve(NonlinearLeastSquaresProblem(err_nls, p0_3b, (trans_KBBaa, params_base, fitdat)), LevenbergMarquardt())
sol3b = gen_sol_pd(opt3b.u, trans_KBBaa, params_base)
prm3b = sol3b.prob.p
@show prm3b.Rp


begin
blankplot_hrC()
modrftplot!(sol3b)
plot!(fitdat)
end

## ------------------
# Tune with model version 2, which is DIC

trans_KBBα = merge(trans_KBB, as((;alpha=TVScale(0.001u"cm^(3//2)") ∘ TVExp())))
# p0_2 = vcat(opt3.u, [0.0])
p0_2 = [3.0, 3.0, -1, -1]
# ub2 = [10.0, 7.0, 7.0, 3.0]
# lb2 = -ub2
# opt2 = solve(OptimizationProblem(objf_2, p0_2, (params_base, fitdat), lb=lb2, ub=ub2), NelderMead())
err_nls_2 = NonlinearFunction{true}(resid_prototype=zeros(num_errs(fitdat))) do resid, fitlog, tpf
    sol = gen_sol_rf_LC2(fitlog, tpf...)
    LyoPronto.err_expT!(resid, sol, tpf[3])
end
@time opt2 = solve(NonlinearLeastSquaresProblem(err_nls_2, p0_2, (trans_KBBα, params_base, fitdat)), LevenbergMarquardt() )
sol2 = gen_sol_rf_LC2(opt2.u, trans_KBBα, params_base)
prm2 = sol2.prob.p

begin
blankplot_hrC(ylim=(-40, 60))
modrftplot!(sol2, trimend=1)
plot!(fitdat)
end
# savefig(plotsdir("M1_LC2.svg"))


# ----------------
# Tune with model version 1
# objf_1 = OptimizationFunction(obj_KBBa1, AutoForwardDiff())
p0_1 = copy(opt2.u)
# @time objf_1(p0_1, (params_base, fitdat))
# ub1 = [10.0, 7.0, 7.0, 3.0]
# lb1 = -ub1
# opt1 = solve(OptimizationProblem(objf_1, p0_1, (params_base, fitdat), lb=lb1, ub=ub1), NelderMead())
# opt1 = solve(OptimizationProblem(objf_1, p0_1, (params_base, fitdat)), NelderMead(), show_trace=true, g_tol=1e-2)

err_nls_1 = NonlinearFunction{true}(resid_prototype=zeros(num_errs(fitdat))) do resid, fitlog, tpf
    sol = gen_sol_rf_LC1(fitlog, tpf...)
    LyoPronto.err_expT!(resid, sol, tpf[3])
end
@time opt1 = solve(NonlinearLeastSquaresProblem(err_nls_1, p0_1, (trans_KBBα, params_base, fitdat)), LevenbergMarquardt() )
sol1 = gen_sol_rf_LC1(opt1.u, trans_KBBα, params_base)
prm1 = sol1.prob.p

begin
blankplot_hrC()
modrftplot!(sol1, trimend=1)
plot!(fitdat)
end

#---------------------
# Fit taken directly from SciRep article
p_sr = @set params_base.Kvwf = 1e-3u"cal/s/K/cm^2"
p_sr = @set p_sr.Bf = 2e7u"Ω/m^2" /0.54 # un-correct efficiency
p_sr = @set p_sr.Bvw = 0.9e7u"Ω/m^2" /0.54
p_sr = @set p_sr.alpha = 0.08u"cm^(3//2)"
p_sr = @set p_sr.eppf = ConstPhysProp(8e-4)
p_sr = @set p_sr.eppvw = 2.4e-2

# prob_sr = ODEProblem(lumped_cap_rf_alt1, LyoPronto.calc_u0(p_sr), (0.0, 1000.0), p_sr)
prob_sr = ODEProblem(lumped_cap_rf_alt1, [5.0, 235.15, 235.15], (0.0, 1000.0), p_sr)
sol_sr = solve(prob_sr, Rosenbrock23(autodiff=AutoFiniteDiff());  callback=end_drying_callback)

# tTf_sr = CSV.read(datadir("sims", "M1srep_Tf.csv"), Table)
# tTvw_sr = CSV.read(datadir("sims", "M1srep_Tvw.csv"), Table)
# t_sr = (tTf_sr.t ∩ tTvw_sr.t)
# Tf_sr = tTf_sr.Tf[tTf_sr.t .∈ [t_sr]]u"°C"
# Tvw_sr = tTvw_sr.Tvw[tTvw_sr.t .∈ [t_sr]]u"°C"
# fit_sr = PrimaryDryFit(t_sr*u"hr", Tf_sr, Tvw_sr, 10.8u"hr")
# plot(fit_sr, trim=1)

# p0_sr = [0, 0.1, 0.1, 0.1]
# transform(trans_KBBα, p0_sr)

# err_nls_sr = NonlinearFunction{true}(resid_prototype=zeros(num_errs(fit_sr))) do resid, fitlog, tpf
#     sol = gen_sol_rf_LC1(fitlog, tpf...)
#     LyoPronto.err_expT!(resid, sol, tpf[3])
# end
# @time opt_sr = solve(NonlinearLeastSquaresProblem(err_nls_sr, p0_sr, (trans_KBBα, p_sr, fit_sr)), LevenbergMarquardt() )
# sol_sr = gen_sol_rf_LC1(opt_sr.u, trans_KBBα, params_base)
# prm_sr = sol_sr.prob.p

begin
blankplot_hrC()
modrftplot!(sol_sr, trimend=1)
plot!(fitdat)
# modrftplot!(sol1, c=:red, trimend=1)
end

# ----------------
# Plot all variants together

# begin
# blankplot_hrC()
# plot!(fitdat, alpha=1.0, c=:black, markers=[:plus :diamond :square :none])
# modrftplot!(sol1,  labsuffix=", LC1",  c=1, markers = :circle)
# modrftplot!(sol1b, labsuffix=", LC1b", c=2, markers=:square)
# modrftplot!(sol2,  labsuffix=", LC2",  c=3, markers=:dtriangle)
# modrftplot!(sol3,  labsuffix=", LC3",  c=4, markers=:utriangle)
# plot!(legend=:outerright)
# end
# savefig(plotsdir("M1_allLC_share.svg"))
# savefig(plotsdir("M1_allLC_share.pdf"))

opts = [opt1, opt2, opt3]
sols = [sol_sr, sol1, sol2, sol3]
errs = sqrt.([obj_expT(sol_sr, fitdat), map(o->sum(abs2, o.resid), opts)...])
prms = [p_sr, prm1, prm2, prm3]
modelnames = ["IC-m", "IC", "DIC", "DIF"]
ticklabs = (1:4, modelnames)

fitdat_less = @set fitdat.Tfs = fitdat.Tfs[2:2]
@reset fitdat_less.t_end = missing
default(:lw, 3)

begin
pl_sr, pl_lc1, pl_lc2, pl_lc3 = map(sols, modelnames) do sol, modname
    pl = blankplot_hrC(ylims=(-42, 40), xlims=(0, 12.45))
    plot!(fitdat_less,nmarks=20, labsuffix=", experiment")
    modrftplot!(sol, trimend=1, sampmarks=false)
    tendplot!(fitdat.t_end, ls=:dash, label="")
    plot!(title=modname, legend=:none, )
end
# plot(pl_lc1, pl_lc1b, pl_lc2, pl_lc3, layout=(2,2)) 
annotate!(pl_sr, 7, -30, Plots.text("end of drying", 10, "Computer Modern"))
plot!(pl_sr, [10, 11.5], [-30, -30], arrow=:arrow, c=:black, linewidth=1, label="")

plot!(pl_sr, left_margin=20Plots.px)
plot!(pl_lc1, left_margin=-10Plots.px, ylabel="", yticks=(-40:20:40, ""))
plot!(pl_lc2, left_margin=-10Plots.px, )
plot!(pl_lc3, left_margin=-10Plots.px, ylabel="", yticks=(-40:20:40, ""))
pl_err = bar(errs, xticks=ticklabs, label="", lw=1, title="RMS Error", yguide=L"$\sqrt{L(\vec{\theta})}$  $[\mathrm{K}]$", left_margin=20Plots.px, ylims=(0,11))
bar!(pl_err, [5], [9.465], xticks=(1:5, vcat(modelnames, ["LS"])))
plot!(pl_err, left_margin=10Plots.px)
pl_labs = deepcopy(pl_sr)
plot!(pl_labs, ylim=(-1,0), xlim=(-1,0), frame=:none, legend=(0, 0.5), ylabel=nothing, xlabel=nothing, title="", legendfontsize=12)
plot(pl_sr, pl_lc1, pl_labs, pl_lc2, pl_lc3, pl_err, bottom_margin=20Plots.px, layout=@layout([a b c{0.3w}; d e f]), size=(800, 600)) 
end
# savefig(plotsdir("M1_allLC_sep.svg"))
# savefig(plotsdir("M1_allLC_sep.pdf"))
savefig(plotsdir("M1_allLC_sep_v3.svg"))

# opts = [opt1, opt1b, opt2, opt3]
# prms = [prm1, prm1b, prm2, prm3]
# ticklabs = (1:4, ["LC1", "LC1b", "LC2", "LC3"])
function row_from_params(prm)
    nt = (a1 = prm.Rp.A1, a2 = prm.Rp.A2)
    for nm in [:alpha, :Bf, :Bvw, :Kvwf]
        nt = merge(nt, (nm => getfield(prm, nm),))
    end
    return nt
end

table = Table(map(row_from_params, prms))
# safesave(plotsdir("M1_tuning_tab.html"), pretty_table(HTML, table))

# table.alpha[4:4] .= 0.0u"cm^1.5"

formatter = (label, unit)-> label *"\n\n"* latexify(unit)
markers = [:circle, :square, :diamond, :hexagon]
colors = [1, 2, 3, 4]
begin
resetfontsizes()
fitsattr = (label="", c=colors, markers=markers, unitformat=latexify, ylabel=" ", left_margin=20Plots.px, widen=1.2, grid=:y, xticks=:none, markersize=5)
pl1 = @df table scatter(:Kvwf .|> u"W/m^2/K",   title=L"K_\mathrm{vw-f}"; fitsattr...)
pl2 = @df table scatter(:Bf,  title=L"B_\mathrm{f}"; yscale=:log10, ylim=(4e6, 2e9), fitsattr...)
hline!([1.5e7], c=:gray, ls=:dash, label="")
# pl2 = @df table scatter(:Bf,  ylabel=L"B_\text{f}"; yticks=8e8:2e8:1.4e9, fitsattr...)
# pl3 = @df table scatter(:Bvw, title=L"B_\mathrm{vw}"; yscale=:log10, ylim=(4e6, 2e7), fitsattr...)
pl3 = @df table scatter(:Bvw, title=L"B_\mathrm{vw}"; yticks=[1e7, 1e8, 1e9], ylim=(4e6,2e9), yscale=:log10, fitsattr...)
hline!([1.2e7], c=:gray, ls=:dash, label="")
plot!(yticks=([1e7, 1e8, 1e9],[]),ylabel="", left_margin=-20Plots.px)
pl4 = @df table scatter(:alpha, c=[1, 1, 1, :black],  title=L"\alpha",  widen=1.2; fitsattr...)
# pl5 = @df table scatter(:a1, c=[:black, 1, :black, :black], ylabel=L"a_1", widen=1.2; fitsattr...)
# pl6 = @df table scatter(:a2, c=[:black, 1, :black, :black], yla3el=L"a_2", widen=1.2; fitsattr...)
for pl in [pl1, pl2, pl3, pl4]
    plot!(pl, xticks=ticklabs, xtick_dir =:none)
end
# plot(pl1, pl4, pl2, pl5, pl3, pl6, layout=(3,2), size=(600,400))
plot(pl1, pl2, pl3, pl4, formatter=x->latexify(x, fmt=FancyNumberFormatter(2, "\\times")), layout=(1,4), size=(800, 200), top_margin=15Plots.px)
end
savefig(plotsdir("M1_allLC_params.svg"))
savefig(plotsdir("M1_allLC_params.pdf"))

# begin
# lr = range(0u"cm", h_f0, length=40) .|> u"cm"
# RpAs = map(prms) do prm
#     A_sub = A_p .+ (ismissing(prm.alpha) ? 0u"cm^2" : prm.alpha*sqrt.(lr)) .|> u"cm^2"
#     prm.Rp.(lr)./A_sub
# end
# rpattr = (seriestype=:samplemarkers, markersize=5, linealpha=0.7)
# pl4 = plot(u"cm", u"hr*Torr/g", xlabel=L"h_\text{d}", ylabel=L"deviation in $R_p / A_\text{sub}$", unitformat=:square)
# for (i, rp) in enumerate(RpAs)
#     plot!(lr, rp.-RpAs[1], label=ticklabs[2][i], offset=i, step=4, m=markers[i]; rpattr...)
# end
# # layout = @layout [[a ;b; c] d{0.7w}]
# # plot(pl1, pl2, pl3, pl4, layout=layout, size=(800,400))
# plot!(legend=:bottomleft)
# plot!(xticks=0:0.3:2.0, size=(400,400))
# end
# savefig(plotsdir("M1_allLC_RpA.svg"))
# savefig(plotsdir("M1_allLC_RpA.pdf"))

# -----------------
# Braatz group analytical model

const AM = MicrowaveLyoModeling.AnalyticalModel

Tm = find_zero( T->LyoPronto.calc_psub(T*u"K")-pch(0), 250)*u"K"
Kv = K_shf_f(pch(0))
Qppp_f = uconvert(u"W/m^3", P_per_vial(0)*prm1.Bf*(2π*LyoPronto.e_0)*f_RF*LyoPronto.eppf(Tm, f_RF))
kf = LyoPronto.k_ice
Tb0 = Tsh.setpts[1]
Tbf = Tsh.setpts[2]
r = Tsh.ramprates[1]

pp = AM.Params(Kv, Qppp_f, kf, h_f0, r, Tb0, Tbf, Tm)
t, s, T = AM.calc_tsT(pp);

begin
pl_all = blankplot_hrC()
@df thm_pd exptfplot!(:t, :T4, nmarks=40, markers=[:circle :square])
@df thm_pd exptvwplot!(:t, :T3, nmarks=40)
tendplot!(fitdat.t_end, ls=:dash, label="")
modrftplot!(sol3, labsuffix=", LC-DIF",trimend=1)
plot!(t, T, c=:green, label=L"$T_f$, TLM")
plot!(legend=:topleft)
plot!([1, 2.7], [-37, -27], arrow=:arrow, c=:gray, linewidth=2, label="")
pl_brtz = plot!(pl_all, u"hr", u"°C"; inset=bbox(0.38, 0.65, 0.3, 0.15), subplot=2)
pl_brtz = pl_all[2]
plot!(pl_brtz, t, T ,label="",  markersize=7, c=:green)
plot!(pl_brtz; ylim=(-40.1, -39.4), xlabel="", ylabel="", )
end
savefig(plotsdir("M1_MIT_compT.svg"))
savefig(plotsdir("M1_MIT_compT.pdf"))
