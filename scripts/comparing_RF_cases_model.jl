using DrWatson
@quickactivate :MicrowaveLyoModeling

using Latexify
resetfontsizes()

plot_defaults_lprf()

casenames = ["M1", "M2", "M3", "M4", "SM"]
casenames_pretty = copy(casenames); casenames_pretty[5] = "SM2"
# casenames_pretty = ["mannitol", "man:suc 2:1", "BSA", "ALV", "mAb:Suc 1:1", "mAb:Suc 1:5", "mAb:Tre 1:1", "mAb:Tre 1:5"]

# Manually added for now.
vialsizes = ["6R", "6R", "6R", "6R", "6R"]
fillvols = [5u"mL", 3u"mL", 3u"mL", 3u"mL", 3u"mL"]

cases = Table((name=casenames, prettyname=casenames_pretty, 
        vialsize=vialsizes, fillvol=fillvols))

fnames_param = [datadir("exp_pro", nm*"_KvRpRF.jld2") for nm in cases.name]
fnames_q = [datadir("exp_pro", nm*"_Q.jld2") for nm in cases.name]

caseparams = map(load, fnames_param)
caseqs = map(load, fnames_q)

function colrename(res)
    rad = get_vial_radii.(cases.vialsize[1])[1]
    Bi = res["Kvwf"]*rad/LyoPronto.k_ice |> NoUnits
    nt = (Kvwf = res["Kvwf"] |> u"W/m^2/K", 
          Bf = res["Bf"], 
          Bvw = res["Bvw"], 
          fBf = res["Bf"]*res["po"].f_RF,
          fBvw = res["Bvw"]*res["po"].f_RF,
          R0 = res["Rp"].R0, 
          A1 = res["Rp"].A1,
          A2 = res["Rp"].A2,
          Kshf = res["Kshf"],
          Tferr = res["Tferr"],
          Tvwerr = res["Tvwerr"],
          terr = res["terr"],
          Bi)
end
allfits = Table(map(colrename, caseparams))


allqs = Table(map(x->(QRFf=x["QRFf"], Qvwf=x["Qvwf"], Qshf=x["Qshf"]), caseqs))
cases_res = Table(cases, allfits, allqs)
@show cases_res.Bi

format_2ln = (label, unit)-> label *"\n\n"* latexify(unit)
format_1ln = (label, unit)-> label *"\n["* latexify(unit) * "]"
ticklabs = (range(1, length(cases)), cases.prettyname)

# begin
# # scatter(Kshfs, label="", ylabel=L"K_\text{sh-f}", unitformat=formatter)
# @df allfits scatter(:Kshf, label="", ylabel="K_{sh-f}", unitformat=latexsquareunitlabel)
# plot!(grid = :y, xtick_dir = :none, tickfontrotation=45, xticks=(1:6, casenames_pretty))
# plot!(xticks = ticklabs, size=(400,400))
# end
# savefig(plotsdir("compare_Kv.svg"))
# savefig(plotsdir("compare_Kv.pdf"))

# begin
# rpattr = (xticks=:none, markersize=6, label="", ywiden=1.2, yscale=:identity, unitformat=format_2ln)
# pl1 = @df allfits scatter(:R0; ylabel=L"R_0", rpattr...)
# pl2 = @df allfits scatter(:A1; ylabel=L"A_1", rpattr...)
# pl3 = @df allfits scatter(:A2; ylabel=L"A_2", rpattr...)
# plot!(grid = :y, xtick_dir = :none, xtickfontrotation=45)
# plot!(xticks=ticklabs, )
# plot(pl1, pl2, pl3, layout=(3,1), link=:x, size=(400,400))
# end
# savefig(plotsdir("compare_Rp_params.svg"))
# savefig(plotsdir("compare_Rp_params.pdf"))

begin
rfattr = (xticks=(1:5, []), grid=:y, c=[1, :black, 1, 1, 1], markersize=5, label="", yscale=:log10, ywiden=1.2, xwiden=1.2, unitformat=format_1ln)
pl1 = @df allfits scatter(:Kvwf, ylabel=L"K_\mathrm{vw-f}",;  rfattr...,c=[1, :black, 1, :black, 1], )
pl2 = @df allfits scatter(:Bf,   ylabel=L"B_\mathrm{f}"; ylims=(4e6, 2e9), rfattr...)
hline!([1.5e7], l=:dash, c=:gray, label="")
pl3 = @df allfits scatter(:Bvw,  ylabel=L"B_\mathrm{vw}"  ; ylims=(4e6,2e9),rfattr...)
hline!([1.2e7], l=:dash, c=:gray, label="")
# pl4 = @df allfits scatter(:f,  ylabel=L"B_\text{vw}"  ; rfattr...)
plot!(xtick_dir =:in, xtickfontrotation=45)
plot!(xticks=ticklabs)
plot!(left_margin=20Plots.px)
pl_RFp = plot(pl1, pl2, pl3, layout=(3,1), link=:x, size=(400,400))
end
# savefig(plotsdir("compare_RFparams.svg"))
# savefig(plotsdir("compare_RFparams.pdf"))

begin
pl1 = @df allfits bar(sqrt.(:Tferr), xticks=:none, lw=1, label="", ylim=(0, 20.5))
plot!(ylabel="RMS error,  \n" * L"$T_\mathrm{f} \quad [\mathrm{K}]$")
pl2 = @df allfits bar(sqrt.(:Tvwerr), xticks=:none, lw=1, label="", c=[1, 2, 1, 1, 1], ylim=(0, 20.5))
plot!(ylabel="RMS error,  \n" * L"$T_\mathrm{vw} \quad [\mathrm{K}]$")
pl3 = @df allfits bar(sqrt.(:terr), xticks=ticklabs, lw=1, label="")
plot!(ylabel="Abs. error,  \n" * L"t_\mathrm{end}\quad [\mathrm{hr}]")
plot!(ylims=(0, 1.8), xgrid=false, xtick_dir=:none)
ple = plot(pl1, pl2, pl3, layout=(3,1), size=(400,400), link=:x)
end
# savefig(plotsdir("compare_RF_err.svg"))
# savefig(plotsdir("compare_RF_err.pdf"))

begin
@df allqs barstackplot(ticklabs[1], :QRFf, :Qvwf, :Qshf, labels=[L"Q_\mathrm{RF-f}" L"Q_\mathrm{vw-f}" L"Q_\mathrm{sh-f}"] )
plot!(xticks = ticklabs)
plot!(yticks = [], ylabel="Relative Magnitude of Heating", left_margin=25Plots.px)
plot!(size=(400, 300), legend=:outerright)
pl_Q = plot!(grid = :y, xtick_dir = :in, tickfontrotation=45)
end
# savefig(plotsdir("compare_Q.svg"))
# savefig(plotsdir("compare_Q.pdf"))

# plot(pl_RFp, pl_Q, layout=(1,2), size=(800, 400))
# savefig(plotsdir("compare_RFQ.svg"))
# savefig(plotsdir("compare_RFQ.pdf"))


plot(pl_RFp, ple, pl_Q, layout=@layout([a b c{0.4w}]) , size=(900, 400), left_margin=28Plots.px)
savefig(plotsdir("compare_RFQe.svg"))
savefig(plotsdir("compare_RFQe.pdf"))
# blankmmplot() = plot(u"mm", u"mm", ylabel="",xlabel="", showaxis=false, grid=false, aspect_ratio=:equal, yscale=:identity)

# begin
# gr()
# figs = map(zip(vialsizes, fillvols, casenames_pretty)) do (v, f, n)
#     voutline, foutline = make_outlines(get_vial_shape(v), f)
#     vshape = Shape(voutline)
#     fshape = Shape(foutline)

#     pl = blankmmplot()
#     plot!(vshape, fillcolor=:lightgray, label="", lw=1)
#     plot!(fshape, fillcolor=:blue, label="", lw=1)
#     plot!(xlims=(-15, 15))
#     plot!(ylims=(0, 50))
#     plot!(xlabel = n, xguidefontrotation=0)
#     plot!(left_margin=-30Plots.px, right_margin=-10Plots.px)
#     return pl
# end
# yvals = [i*u"cm" for i in 0:1:5]
# ylabs = ["$i cm" for i in 0:1:5]
# # plot!(figs[1], showaxis=:y, yticks=(yvals, ylabs), left_margin=0Plots.px)
# plot!(figs[6], showaxis=:y, yticks=(yvals, ylabs), left_margin=0Plots.px)
# plot(figs..., link=:x, layout=grid(1,length(casenames)), size=(600, 200), top_margin=0Plots.px, bottom_margin=30Plots.px)
# end
# savefig(plotsdir("compare_vialfill.svg"))
# savefig(plotsdir("compare_vialfill.pdf"))
