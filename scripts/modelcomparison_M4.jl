plot_defaults_mlm()

M4 = load(datadir("exp_pro", "M4_processed.jld2"))
@unpack thm_pd, lyo_pd, fitdat, P_per_vial = M4

tstops = P_per_vial.t

po_M1 = load(datadir("exp_pro", "M1_KvRpRF.jld2"))["po"]


po_M4 = @set po_M1.P_per_vial = P_per_vial
@reset po_M4.Tsh = RampedVariable(u"K".([-40u"°C", 10u"°C"]), 0.5u"K/minute")
@reset po_M4.pch = RampedVariable(lyo_pd.pch_sp[1])
Vfill = 3u"mL"
@reset po_M4.hf0 = Vfill/po_M1.Ap
@reset po_M4.mf0 = Vfill*po_M1.ρsolution
@reset po_M4.f_RF = 18u"GHz"

tsol = solve(ODEProblem(po_M4), Rodas3(), tstops=ustrip.(u"hr", tstops))
modrftplot(tsol)
plot!(fitdat)

# Not adequate parameters. At different field strength, so this makes sense, probably
using LyoPronto.TransformVariables: logit
trans = as((Bf = TVScale(7e7u"Ω/m^2") ∘ TVScale(1e2) ∘ TVLogistic() ∘ TVShift(logit(1e-2)),
            Bvw = TVScale(1e7u"Ω/m^2") ∘ TVScale(1e2) ∘ TVLogistic() ∘ TVShift(logit(1e-2)) ))
p0 = [1.0, -1.0]
# objf = OptimizationFunction(obj_pd, AutoForwardDiff(chunksize=2))
# @time objf(p0, (trans, po_M4, fitdat))
nls_M4 = NonlinearFunction{true, SciMLBase.FullSpecialize}(nls_pd!, resid_prototype=zeros(num_errs(fitdat)))
opt = solve(NonlinearLeastSquaresProblem(nls_M4, p0, (trans, po_M4, fitdat)), LevenbergMarquardt())
prof_RF = gen_sol_pd(opt.u, trans, po_M4)

begin
plT = blankplot_hrC(legend=:topleft)
@df thm_pd exptfplot!(:t, :T3, sampmarks=true, linealpha=0.5, nmarks=30)
@df thm_pd exptvwplot!(:t, :T1, :T4, sampmarks=true, linealpha=0.5, nmarks=40)
modrftplot!(prof_RF, markeralpha=0, evensample=false)
plot!(po_M4.Tsh, tmax=fitdat.t_end, c=:black, label=L"T_\mathrm{sh}")
tendplot!(fitdat.t_end, ls=:dash, label="")
annotate!(3.8, -30, Plots.text("end of drying", 12, "Computer Modern"))
plot!([fitdat.t_end-1u"hr", fitdat.t_end-0.2u"hr"], [-30, -30], arrow=:arrow, c=:black, linewidth=1, label="")
# end
# savefig(plotsdir("M4_compT.svg"))
# savefig(plotsdir("M4_compT.pdf"))

tt = range(0u"hr", 7u"hr", length=3600)
pp = [P_per_vial(ti) for ti in tt]
# begin
plp = plot(u"hr", u"W", xlabel="Time", ylims=(0, 0.7))
tendplot!(fitdat.t_end, ls=:dash, label="")
# for (ti, to) in zip(tstops[1:2:end], tstops[2:2:end])
#     vspan!([ti, to], lw=0, c=:red, alpha=0.4, label="")
# end
plot!(tt, pp, fillrange=0, lw=0, c=:red, label="")
plot!(; yticks=[0.0, 0.7], widen=false, ylabel="RF Power\n[W/vial]", yunitformat=(l,u)->l, top_margin=-30Plots.px)
# plot!(yticks=[], ylabel="RF", size=(600, 100), xtickdir=:out, widen=false)
plot!(plT, xlabel="", xwiden=false, legend=:bottomright, xformatter=x->"", left_margin=30Plots.px, )
plTp = plot(plT, plp, link = :x, xlims=(0,8.5), layout=@layout([a; b{0.1h}]), size=(600,400), bottom_margin=20Plots.px)
# savefig(plotsdir("M4_power.svg"))
# savefig(plotsdir("M4_power.pdf"))

plq = qplotrf(prof_RF, ordering=[3,2,1])
plot!(widen=false, ylim=(0, 0.6), left_margin=20Plots.px)
plot(plTp, plq, layout=@layout([a  b{0.35w}]), size=(1000, 400), bottom_margin=20Plots.px)
end
savefig(plotsdir("M4_T_q_combine.svg"))
savefig(plotsdir("M4_T_q_combine.pdf"))

# Save fit results to a file
save_fitresults(opt, "M4")