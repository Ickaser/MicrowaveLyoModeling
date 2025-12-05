
plot_defaults_mlm()

# -----------------------

begin
simgridsize = (41, 31)
base_props = LevelSetSublimation.base_props

# Get some stuff from LyoProntoNIIMBLRF
lcfitparams = load(datadir("exp_pro", "M1_KvRpRF.jld2"))
@unpack Kvwf, Bf, Bvw = lcfitparams
expdat = load(datadir("exp_pro", "M1_processed.jld2"))
thmdat = expdat["thm_pd"]
t_end = expdat["t_end"]

# ---- Properties which do not change in time
vialsize = "6R"
fillvol = 5u"mL"
# Mass transfer
c_solid = 0.05u"g/mL"
ρ_solution = 1.0u"g/mL"
ϵ = (ρ_solution-c_solid)/ρ_solution 
κ = 0.0u"m^2" # no viscous-regime flow allowed
# R_p values to mass transfer
Rp0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16u"cm*hr*Torr/g"
Tguess = 260u"K"
l_base = sqrt(u"R"*Tguess/base_props.Mw) / A1
l = fill(l_base, simgridsize)
# l .*= 2.0 .^range(-1, 1, length=simgridsize[1]) # Big pores at the outer radius
l .*= range(0.25, 1.5, length=simgridsize[1]) # Big pores at the outer radius
# Heat transfer
kd = LyoPronto.k_sucrose * (1-ϵ)
m_v = LyoPronto.get_vial_mass(vialsize)
A_v = π*LyoPronto.get_vial_radii(vialsize)[2]^2
# Microwave
B_d = 0.0u"Ω/m^2"
tcprops = TimeConstantProperties(ϵ, l, κ, Rp0, kd, Kvwf, m_v, A_v, B_d, Bf, Bvw)

# ------- Properties which may change in time
f_RF = RampedVariable(8.0u"GHz")
pch = RampedVariable(100.0u"mTorr")
Tsh = RampedVariable(uconvert.(u"K", [-40.0, 10]*u"°C"), 1u"K/minute")
# P_per_vial = RampedVariable(10u"W"/17 * 0.54)
P_per_vial = RampedVariable(0.0u"W")
# Heat transfer coefficient as function of pressure
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)
tvprops = TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)

paramsd = base_props, tcprops, tvprops

# ------------
# Assemble simulation

# config = Dict{Symbol, Any}()
# @pack! config = paramsd, vialsize, fillvol, simgridsize
config = @ntuple paramsd vialsize fillvol simgridsize
config = merge(config, (time_integ=Val(:dae_then_exp),))
end

# Run simulation
@info "before solve"

# @time res, fname = produce_or_load(sim_from_dict, config; filename=hash, verbose=true, tag=true, prefix=datadir("sims", "M1_porevar"))

@time res = sim_from_dict(config; verbose=true)
fname = datadir("sims", "M1_porevar_"*string(hash(config), base=10)*".jld2")
safesave(fname, res)

@info "after solve"

# fname = datadir("sims", "M1_"*string(hash(config), base=10)*".jld2")
# res = load(fname)

sim = res["sim"]
@unpack dom = sim

# --------------

# Plot things of interest

# ---------------------------
# Plots

# summaryT(sim, tstart=0.01, tend=0.99, layout=(2,3), clims=(-30, 50))
# plot!(size=(800,300), bottom_margin=10Plots.px)
# savefig(plotsdir("M1_LSS_summary.svg"))
# savefig(plotsdir("M1_LSS_summary.pdf"))

locs = [(0.0, 0.02), (0.6, 0.05)]#, (0.2, 0.5)]
vtmarks = [:circle, :utriangle, :circle]
# vtmarks = [:diamond, :utriangle, :square, :circle]

resetfontsizes(); scalefontsizes(1.2)
# Plot pore size variation
begin
hextend = 1.10
vextend = 1.05
walls = Plots.Shape([(-dom.rmax*hextend, dom.zmax+dom.dz*2), 
            (-dom.rmax*hextend, -(vextend-1)*dom.zmax),
            ( dom.rmax*hextend, -(vextend-1)*dom.zmax),
            ( dom.rmax*hextend, dom.zmax+dom.dz*2),
            ( dom.rmax+dom.dr/2, dom.zmax+dom.dz*2),
            ( dom.rmax+dom.dr/2, -dom.dz/2),
            (-dom.rmax-dom.dr/2, -dom.dz/2),
            (-dom.rmax-dom.dr/2, dom.zmax+dom.dz*2)])
ll = l .|> u"μm"
pl_l = plot(aspect_ratio=:equal, showaxis=false, grid=false, widen=false, xlabel="Imposed pore variation")
heatmap!(dom.rgrid .- dom.rmax, dom.zgrid, ll[end:-1:begin, :]', cmap=:batlow)
heatmap!(dom.rgrid, dom.zgrid, ll', cmap=:batlow, colorbar_title="Pore size\n")
plot!(walls, c=:gray, lw=1, label="")
plot!(bottom_margin=-10Plots.px, right_margin=0Plots.px, left_margin=0Plots.px, top_margin=-10Plots.px)
# plot!([-dom.rmax, -dom.rmax, dom.rmax, dom.rmax], [dom.zmax, -dom.dz/2, -dom.dz/2, dom.zmax], c=:black, lw=2, label="")
# end

# begin
tf = sim.sol.t[end]
tsums = range(5, 25, length=3)
pl_sum = summaryT(sim, tstart=first(tsums)*3600/tf, tend=last(tsums)*3600/tf, layout=(3,1))
for i in 1:3
    # placethermocouples!(pl_sum[i], dom, locs, c=palette(:Oranges_4)[4:-1:2], markers=vtmarks, label="", markersize=8);
    # placethermocouples!(pl_sum[i], dom, [(0.95, 0.9)], msc=palette(:Oranges_4)[4:4], c=:white, markers=[:circle], label="", markersize=8, msw=4);
    plot!(pl_sum[i], LevelSetSublimation.PlaceThermocouples((dom, locs)), c=palette(:Oranges_4)[4:-1:2], markers=vtmarks, label="", markersize=8);
    plot!(pl_sum[i], LevelSetSublimation.PlaceThermocouples((dom, [(0.95, 0.9)])), msc=palette(:Oranges_4)[4:4], c=:white, markers=[:circle], label="", markersize=8, msw=4);
end
plot!(pl_sum[3], cbar_title="\nTemperature [°C]", left_margin=20Plots.px, right_margin=40Plots.px)
# labs = vcat([L"$T$, model "*i for i in ["bottom","corner","center"]], L"$T_\mathrm{vw}$, model")
labs = vcat([L"$T$, model "*i for i in ["bottom","offset"]], L"$T_\mathrm{vw}$, model")
pl_T = blankplothrC( )
plot!(Tsh, c=:black, tmax=28u"hr", label=L"T_\mathrm{sh}")
# @df thmdat exptfplot!(:t, :T4, nmarks=20)
# @df thmdat exptvwplot!(:t, :T3, nmarks=20, msw=3, thickness_scaling=1.0)
vt_plot!(sim, locs; labels=permutedims(labs), markers=permutedims(vtmarks), step=40, samplemarkers=true, linealpha=0.8)
for t in tsums
    plot!([t, t], [-30, -42], label="", c=:black, lw=1, arrows=true)
end
# tendplot!(t_end, label="")
plot!(legend=:top, bottom_margin=15Plots.px, left_margin=15Plots.px)
# savefig(plotsdir("porevar_multiT.svg"))
# savefig(plotsdir("porevar_multiT.pdf"))
# plot!(pl_sum, margin_right=26Plots.px)
plot(pl_T, pl_l, pl_sum, layout=@layout([[a; b{0.3h}] c]), size=(700,500))
end
savefig(plotsdir("porevar_allLSS.svg"))
savefig(plotsdir("porevar_allLSS.pdf"))

# savefig(plotsdir("M1_multiT_small.svg"))
# savefig(plotsdir("M1_multiT_small.pdf"))


begin
pl_vtloc, T = plotframe(10*60*60, sim)
placethermocouples!(sim.dom, locs, c=palette(:Oranges_4)[4:-1:2], markers=vtmarks, label="", markersize=8);
# texts = [text(L"$T_\textrm{f%$i}$", 10) for i in 1:3]
# labelthermocouples!(sim.dom, locs, texts)
plot!(xlabel="", cbar=nothing)
# annotate!(pl, dom.rmax.*locx, dom.zmax.*locy, texts)
plot!(size=(240,200), colorbar_title="T [°C]")
end
# savefig(pl, plotsdir("M1_virtual_thermocouple_small.svg"))
# savefig(pl, plotsdir("M1_virtual_thermocouple_small.pdf"))

# --------- Compute effective Rp

t_Rp, hd_eff, Rp_eff, Asub = get_eff_Rp(sim);
Rp_orig = @. Rp0 + A1*hd_eff
relerr = (Rp_eff .- Rp_orig)./Rp_orig
begin
pl1 = plot(u"cm", u"cm^2*Torr*hr/g", ylabel="R_p", unitformat=latexify)
plot!(hd_eff, Rp_eff, seriestype=:samplemarkers, marker=:circle, label=L"LS, effective $R_p$ with $l(r,z)$", ylims=(0, 40))
plot!(hd_eff, Rp_orig, c=:black, label=L"Original $R_p$")
plot!(xlabel=nothing, xticks=(0:0.5:1.5, []), bottom_margin=-20Plots.px)
pl2 = plot(u"cm", NoUnits, ylabel="Relative\nDiff.", xlabel=LaTeXString("h_d"), unitformat=latexify)
hline!([0], label="", c=:black, lw=1)
plot!(xticks=(0:0.5:1.5, []), bottom_margin=5Plots.px)
plot!(hd_eff, relerr; c=1, ylim=(-0.18, 0.10), label="",)
pl3 = plot(u"cm", u"cm^2", ylabel="A_\\mathrm{sub}", xlabel=LaTeXString("h_d"), unitformat=latexify)
plot!(hd_eff, Asub, ylim=(3.1, 3.35), label="")
plot!(xticks=(0:0.5:1.5, []))
plot(pl1, pl2, pl3, layout=@layout([a{0.5h}; b{0.2h}; c]), link=:x, xlim=extrema(hd_eff), size=(600, 600))
end
savefig(plotsdir("porevar_curvatureRp.svg"))
savefig(plotsdir("porevar_curvatureRp.pdf"))

# -------- LCLyo model

A_p, A_v = π.* get_vial_radii("6R").^2
c_solid = 0.05u"g/mL"
ρ_solution = 1.0u"g/mL"
LCparams = ParamObjRF((
    (RpFormFit(Rp0, A1, 0u"1/cm"), fillvol/A_p, c_solid, ρ_solution),
    (Kshf, A_v, A_p),
    (pch, Tsh, P_per_vial),
    (fillvol*ρ_solution, base_props.Cpf, m_v, base_props.cp_vw),
    (f_RF(0), base_props.εppf, base_props.εpp_vw),
    (Kvwf, Bf, Bvw)
))

lcprob = ODEProblem(LCParams)
lclyo = solve(lcprob, Rodas3())

begin
blankplothrC()
@df thmdat exptfplot!(:t, :T1, :T4, marker=:circle, )
@df thmdat exptvwplot!(:t, :T3, marker=:square, )
modrftplot!(lclyo, labels=permutedims([i*", LC" for i in ["\$T_f\$", "\$T_{vw}\$"]]))
plot!(size=(400, 200), legend=:outerright)
end
# savefig(plotsdir("M1_lc_small.svg"))
# savefig(plotsdir("M1_lc_small.pdf"))

begin
pl_comp = blankplothrC()
@df thmdat exptfplot!(:t, :T1, :T4, marker=:circle, )
@df thmdat exptvwplot!(:t, :T3, marker=:square, )
vt_plot!(sim, locs; labels=permutedims(labs), markers=permutedims(vtmarks), step=40)
modrftplot!(lclyo, c=:green, linealpha=0.5, labels=permutedims([i*", LC" for i in ["\$T_f\$", "\$T_{vw}\$"]]))
tendplot!(t_end)
plot!(legend=:outerright, size=(600, 200), bottom_margin=15Plots.px, left_margin=15Plots.px)
plot(pl_comp, pl_vtloc, layout=@layout([a{0.8w} b]), size=(800,300))
end