# NOTE: Although this generates a result, it is not very rigorous.
# I have to assume things like Rp for a trehalose-mAb solution,
# as well as chamber pressure when only Pirani is measured.
# What it does do is demonstrate that figure 4 from Park et al. 2021 is way off.

using DataInterpolations
using Roots
const LP = LyoPronto
using Latexify
using NonlinearSolve

plot_defaults_mlm()
set_default(labelformat=:square)
default(:yunitformat, :square)

# ------------------------------------------------
# Figure 1a from Gitter 2019
MWdat_a = map(row->(t=row.t*u"hr", P=row.P*u"W"), CSV.read(datadir("exp_raw", "gitter2019_fig1a", "gitter2019_MW.csv"), Table))
Tdat_a = map(row->(t=row.t*u"hr", T=row.T*u"°C"), CSV.read(datadir("exp_raw", "gitter2019_fig1a", "gitter2019_T.csv"), Table))
pdat_a = map(row->(t=row.t*u"hr", pch=row.pch*u"μbar"|>u"mTorr"), CSV.read(datadir("exp_raw", "gitter2019_fig1a", "gitter2019_pch.csv"), Table))

P_a = LinearInterpolation(MWdat_a.P, MWdat_a.t, extrapolation_right=ExtrapolationType.Constant)
DataInterpolations.integral(P_a, P_a.t[end])/P_a.t[end] # Average power: 39W
pch_a = LinearInterpolation(pdat_a.pch, pdat_a.t, extrapolation_right=ExtrapolationType.Constant)
# plot(P_a, ylim=(0, 120))
# plot(pch_a)

fitdat_a = @df Tdat_a PrimaryDryFit(:t[:T .<-20u"°C"], :T[:T .< -20u"°C"] .|> u"K", Tvws=maximum(Tdat_a.T)|>u"K")
fitdat_a_time = @df Tdat_a PrimaryDryFit(:t[:T .<-20u"°C"], :T[:T .< -20u"°C"] .|> u"K", maximum(Tdat_a.T)|>u"K", (6.0u"hr", 10.0u"hr"))
# plot(fitdat_a_time)

# ----------------------------------------
# # Figure 1b from Gitter 2019
# Load data into memory
MWdat_b = map(row->(t=row.t*u"hr", P=row.P*u"W"), CSV.read(datadir("exp_raw", "gitter2019_fig1b", "gitter2019_MW1b.csv"), Table))
Tdat_b = map(row->(t=row.t*u"hr", T=row.T*u"°C"), CSV.read(datadir("exp_raw", "gitter2019_fig1b", "gitter2019_T1b.csv"), Table))
pdat_b = map(row->(t=row.t*u"hr", pch=row.pch*u"μbar"|>u"mTorr"), CSV.read(datadir("exp_raw", "gitter2019_fig1b", "gitter2019_pch1b.csv"), Table))

fitdat_b = PrimaryDryFit(Tdat_b.t, Tdat_b.T .|> u"K", Tvws=maximum(Tdat_b.T)|>u"K") 
fitdat_b_time = PrimaryDryFit(Tdat_b.t, Tdat_b.T .|> u"K", maximum(Tdat_b.T)|>u"K", (9.0u"hr", 13.5u"hr")) # Bracket the  end time
plot(fitdat_b_time, nmarks=40)

# pch = RampedVariable(minimum(pdat.pch))
pch_b = LinearInterpolation(pdat_b.pch, pdat_b.t, extrapolation_right=ExtrapolationType.Constant)
plot(pch_b)
# @df pdat plot!(:t, :pch, ylim=(0, NaN))
P_b = LinearInterpolation(MWdat_b.P , MWdat_b.t; extrapolation_right=ExtrapolationType.Constant) # Unknown number of vials
DataInterpolations.integral(P_b, P_b.t[end])/P_b.t[end] # Average power: 58W

# -----------------------------------
# Figure 1c from Gitter 2019
MWdat_c = map(row->(t=row.t*u"hr", P=row.P*u"W"), CSV.read(datadir("exp_raw", "gitter2019_fig1c", "MW.csv"), Table))
sort!(MWdat_c, by=r->r.t)
T1dat_c = map(row->(t=row.t*u"hr", T=row.T*u"°C"), CSV.read(datadir("exp_raw", "gitter2019_fig1c", "TC1.csv"), Table))
T2dat_c = map(row->(t=row.t*u"hr", T=row.T*u"°C"), CSV.read(datadir("exp_raw", "gitter2019_fig1c", "TC2.csv"), Table))
pdat_c = map(row->(t=row.t*u"hr", pch=row.pch*u"μbar"|>u"mTorr"), CSV.read(datadir("exp_raw", "gitter2019_fig1c", "pch.csv"), Table))

P_c = LinearInterpolation(MWdat_c.P, MWdat_c.t, extrapolation_right=ExtrapolationType.Constant)
DataInterpolations.integral(P_c, P_c.t[end])/P_c.t[end] # Average power: 83W
pch_c = LinearInterpolation(pdat_c.pch, pdat_c.t, extrapolation_right=ExtrapolationType.Constant)
# plot(P_c, ylim=(0, 120))
# plot(pch_c)

fitdat_c = PrimaryDryFit(T1dat_c.t, (T2dat_c.T[T1dat_c.t.<4u"hr"].|>u"K", ), Tvws=maximum(T1dat_c.T)|>u"K")
fitdat_c_time = PrimaryDryFit(T1dat_c.t, (T2dat_c.T[T1dat_c.t.<4u"hr"].|>u"K", ), maximum(T1dat_c.T)|>u"K", (4.5u"hr", 6.0u"hr"))
# plot(fitdat_c_time)

# Get M1 params, for convenience
M1 = load(datadir("exp_pro", "M1_KvRpRf.jld2"))

begin
# --------------------------
# # Set up parameters: all non-fit parameters are in common across experiments, conveniently.
# Vial parameters
vialsize = "10R"
rad_i, rad_o = get_vial_radii(vialsize)
A_p = π*rad_i^2  # cross-sectional area inside the vial
A_v = π*rad_o^2 # vial bottom area
m_v = get_vial_mass(vialsize)
# Formulation and fill
c_solid = 0.15u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 0.8u"cm^2*hr*Torr/g" # Total guess
A1 = 30.0u"cm*hr*Torr/g" # Total guess
A2 = 1.0u"1/cm" # Total guess
Rp = RpFormFit(R0, A1, A2)
Vfill = 2.3u"mL"
# Heat transfer: ambient chamber wall with some radiation
K_shf_f = ConstPhysProp(1u"W/m^2/K") # Total guess
Tsh = RampedVariable(283.15u"K") # Assumed ambient chamber wall
# Geometry
h_f0 = Vfill/A_p
m_f0 = Vfill * ρ_solution
# RF fit parameters, taken from M1
@unpack Bf, Bvw, Kvwf = M1
# Controllable inputs
f_RF = 2.45u"GHz"

params_a = ParamObjRF((
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf_f, A_v, A_p),
    (pch_a, Tsh, P_a),
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl),
    (f_RF, LyoPronto.eppf, LyoPronto.epp_gl),
    (Kvwf, Bf, Bvw),
))

params_b = @set params_a.P_per_vial = P_b
@reset params_b.pch = pch_b
params_c = @set params_a.P_per_vial = P_c
@reset params_c.pch = pch_c
end

# ---------------- Fit setup
trans_KBB = as((Kvwf=TVScale(Kvwf) ∘ TVExp(),
                Bf = TVScale(Bf) ∘ TVShift(0.01) ∘ TVExp(),
                Bvw = TVScale(Bvw) ∘ TVExp(),))
trans_Rp = Rp_transform_basic(R0, A1, A2)
trans_Kshf = ConstWrapTV() ∘ TVScale(10.0u"W/m^2/K") ∘ TVShift(0.001) ∘ TVLogistic()
trans_KKBBRp = as(merge(trans_KBB.transformations, (;Kshf=trans_Kshf), trans_Rp.transformations))

gensol = (x,tpfu)->gen_sol_pd(x, tpfu[1:3]...; u0=tpfu[4])

u0_a = ustrip.([u"g", u"K", u"K"], [m_f0, fitdat_a.Tfs[1][1], fitdat_a.Tfs[1][1]])
u0_b = ustrip.([u"g", u"K", u"K"], [m_f0, fitdat_b.Tfs[1][1], fitdat_b.Tfs[1][1]])
u0_c = ustrip.([u"g", u"K", u"K"], [m_f0, fitdat_c.Tfs[1][1], fitdat_c.Tfs[1][1]])
p0 = [0.5, -5, -5.5]
p0_KRp = vcat(p0, [-0.1, 0.0, 0.0, 1.0])

nls_g = (du,x,tpfu)->LyoPronto.err_expT!(du, gensol(x, tpfu), tpfu[3])
err_nls_a = NonlinearFunction{true}(nls_g, resid_prototype=zeros(num_errs(fitdat_a)))
err_nls_a_time = NonlinearFunction{true}(nls_g, resid_prototype=zeros(num_errs(fitdat_a_time)))
err_nls_b = NonlinearFunction{true}(nls_g, resid_prototype=zeros(num_errs(fitdat_b)))
err_nls_b_time = NonlinearFunction{true}(nls_g, resid_prototype=zeros(num_errs(fitdat_b_time)))
err_nls_c = NonlinearFunction{true}(nls_g, resid_prototype=zeros(num_errs(fitdat_c)))
err_nls_c_time = NonlinearFunction{true}(nls_g, resid_prototype=zeros(num_errs(fitdat_c_time)))

# ------------ Fit to figure 1a
opt_a_rf = solve(NonlinearLeastSquaresProblem(err_nls_a, p0, (trans_KBB, params_a, fitdat_a, u0_a)), LevenbergMarquardt())
sol_a = gensol(opt_a_rf.u, (trans_KBB, params_a, fitdat_a, u0_a))
prm_a = sol_a.prob.p
sol_a.prob

begin
pl = blankplot_hrC()
plot!(fitdat_a, labsuffix=", experiment")
modrftplot!(sol_a, trimend=0, sampmarks=false)
end

# opt_b = solve(NonlinearLeastSquaresProblem(err_nls, p0_Rp, (trans_KBBRp, params_base, fitdat)), LevenbergMarquardt())
p0_KRp = vcat([-0.5, -2, -4], [-2.1, -1.0, 0.0, 1.0])
p0_KRp = vcat([-0.5, -2, -4], [-1.1, -1.0, 0.0, 1.0])
tsol = gensol(p0_KRp, (trans_KKBBRp, params_a, fitdat_a_time, u0_a))
# obj_expT(tsol, fitdat_a_time, verbose=true)
qplotrf(tsol)
begin
modrftplot(tsol)
plot!(fitdat_a_time)
end

opt_a = solve(NonlinearLeastSquaresProblem(err_nls_a_time, p0_KRp, (trans_KKBBRp, params_a, fitdat_a_time, u0_a)), LevenbergMarquardt())
sol_a = gensol(opt_a.u, (trans_KKBBRp, params_a, fitdat_a_time, u0_a))
prm_a = sol_a.prob.p
# @info "a fit" transform(trans_KKBBRp, opt_a.u)
# obj_expT(tsol, fitdat_a_time, verbose=true)
# obj_expT(sol_a, fitdat_a_time, verbose=true)

qplotrf(sol_a)
# ----------
# Investigate

# Qcontrib = map(sol_a.t) do ti
#     LyoPronto.lumped_cap_rf!(fill(0.0, 3), sol_a(ti), prm_a, ti, Val(true))
# end
# Qcontrib = hcat(Qcontrib...)
# Qsub = Qcontrib[1,:]
# Qshf = Qcontrib[2,:]
# Qvwf = Qcontrib[3,:]
# QRFf = Qcontrib[4,:]
# QRFvw = Qcontrib[5,:]

begin
pl = blankplot_hrC()
plot!(fitdat_a_time, labsuffix=", experiment")
modrftplot!(sol_a, sampmarks=false)
end

begin
pl_T = blankplot_hrC(xwiden=false)
@df Tdat_a exptfplot!(:t, :T, nmarks=40, label=L"$T_\textrm{f}$, Gitter et al. 2019")
modrftplot!(sol_a, labsuffix=", LC-DIF",)
plot!(legend=:topleft, )
# plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)
plot!(xlims=(0, 10.4))
plot!(xticks=0:2:10, bottom_margin=20Plots.px)
pl_q = qplotrf(sol_a, ordering=[1,2,3])
plot!(ylabel="Vial\nHeating\n", ylims=(0, 0.31), xlim=(0, NaN), legend=(0.9, 1))
plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)
pl_mw = plot(xunit=u"hr", xwiden=false, ylabel="MW\nPower\n")
plot!(MWdat_a.t, MWdat_a.P, c=:red, label="")
plot!(ylim=(0, 55), yticks=0:15:60)
plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)
pl_p = plot(xlabel="Time", xunit=u"hr", xwiden=false, ylabel=L"$p_{ch}$"*"\n")
plot!(pdat_a.t, pdat_a.pch, c=3, label="")
plot!(xticks=0:2:10, bottom_margin=20Plots.px)
plot!(ylims=(0, 22), unitformat=:square)
pl_1a = plot(pl_T, pl_q, pl_mw, pl_p; size=(900,400), layout=@layout([a [b; c; d]]), link=:x, left_margin=20Plots.px)
end
savefig(plotsdir("gitter2019_1a.svg"))
savefig(plotsdir("gitter2019_1a.pdf"))


# ------------ Fit to figure 1b
opt_b_rf = solve(NonlinearLeastSquaresProblem(err_nls_b, p0, (trans_KBB, params_b, fitdat_b, u0_b)), LevenbergMarquardt())
sol_b = gensol(opt_b_rf.u, (trans_KBB, params_b, fitdat_b, u0_b))
prm_b = sol_b.prob.p
sol_b.prob

begin
pl = blankplot_hrC()
plot!(fitdat_b, labsuffix=", experiment")
modrftplot!(sol_b, trimend=0, sampmarks=false)
end

# opt_b = solve(NonlinearLeastSquaresProblem(err_nls, p0_Rp, (trans_KBBRp, params_base, fitdat)), LevenbergMarquardt())
p0_KRp = vcat(opt_b_rf.u, [0.1, 1.0, 2.0, 2.0])
opt_b = solve(NonlinearLeastSquaresProblem(err_nls_b_time, p0_KRp, (trans_KKBBRp, params_b, fitdat_b_time, u0_b)), LevenbergMarquardt())
sol_b = gensol(opt_b.u, (trans_KKBBRp, params_b, fitdat_b_time, u0_b))
prm_b = sol_b.prob.p
# @info "b fit" transform(trans_KKBBRp, opt_b.u)

qplotrf(sol_b)

begin
pl = blankplot_hrC()
plot!(fitdat_b_time, labsuffix=", experiment")
modrftplot!(sol_b, sampmarks=false)
end

begin
pl_T = blankplot_hrC(xwiden=false)
@df Tdat_b exptfplot!(:t, :T, nmarks=40, label=L"$T_\textrm{f}$, Gitter et al. 2019")
modrftplot!(sol_b, labsuffix=", LC-DIF",)
plot!(legend=:topleft, )
plot!(xlims=(0, 11.9), xticks=0:2:10, bottom_margin=20Plots.px)
# plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)
pl_q = qplotrf(sol_b, ordering=[1,2,3])
plot!(ylabel="Vial\nHeating\n", ylims=(0, 0.31), xlim=(0, NaN), legend=:topright)
plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)
pl_mw = plot(xunit=u"hr", xwiden=false, ylabel="MW\nPower\n")
plot!(MWdat_b.t, MWdat_b.P, c=:red, fillrange=[0], lw=0, label="")
plot!(ylim=(0, 105), )
plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)
pl_p = plot(xlabel="Time", xunit=u"hr", xwiden=false, ylabel=L"$p_{ch}$"*"\n")
plot!(pdat_b.t, pdat_b.pch, c=3, label="")
plot!(ylims=(0, 17), unitformat=:square)
plot!(xlims=(0, 11.9), xticks=0:2:10, bottom_margin=20Plots.px)
pl_1b = plot(pl_T, pl_q, pl_mw, pl_p; size=(800,400), layout=@layout([a [b; c; d]]), link=:x, left_margin=20Plots.px)
end
savefig(plotsdir("gitter2019_1b.svg"))
savefig(plotsdir("gitter2019_1b.pdf"))

# ------------ Fit to figure 1c

opt_c_rf = solve(NonlinearLeastSquaresProblem(err_nls_c, p0, (trans_KBB, params_c, fitdat_c, u0_c)), LevenbergMarquardt())
sol_c = gensol(opt_c_rf.u, (trans_KBB, params_c, fitdat_c, u0_c))

begin
pl = blankplot_hrC()
plot!(fitdat_a, labsuffix=", experiment")
modrftplot!(sol_a, trimend=0, sampmarks=false)
end

p0_KRp = vcat(opt_c_rf.u, [0.1, 1.0, 2.0, 2.0])
opt_c = solve(NonlinearLeastSquaresProblem(err_nls_c_time, p0_KRp, (trans_KKBBRp, params_c, fitdat_c_time, u0_c)), LevenbergMarquardt())
sol_c = gensol(opt_c.u, (trans_KKBBRp, params_c, fitdat_c_time, u0_c))
prm_c = sol_c.prob.p
# @info "c fit" transform(trans_KKBBRp, opt_c.u)

qplotrf(sol_c)

begin
pl = blankplot_hrC()
plot!(fitdat_c_time, labsuffix=", experiment")
modrftplot!(sol_c, sampmarks=false)
end

# ----------------- Figure 1c, Braatz group comparison
t3 = [0, 1, 3.9]u"hr" # Manually set to match figure
T3 = [236.85, 256.15, 256.15]u"K" # Manually set to match figure

begin
pl_T = blankplot_hrC(xwiden=false)
# plot!(Tsh, tmax=8u"hr", c=:black, label=L"T_\textrm{sh}")
exptfplot!(T1dat_c.t, T1dat_c.T, T2dat_c.T, nmarks=30, label=[L"$T_\textrm{f}$, Gitter et" "    al., 2019"])
modrftplot!(sol_c, labsuffix=", LC-DIF")
# modrftplot!(sol_c, label=["f, LC-DIF" "vw, LC-DIF"])
plot!(xticks=(0:2:10), bottom_margin=20Plots.px)

# plot!(t1, T1, c=:green, label=L"$T_\textrm{f}$, TLM with $T_{sub}$ from $p_{ch}$")
# plot!(t2, T2, c=:green, label="")
# plot!(t2, T2, c=:green, label=L"$T_\textrm{f}$, TLM ")
plot!(t3, T3, c=:purple, label=L"$T_\textrm{f}$, TLM")
plot!(legend=:topleft, )

pl_q = qplotrf(sol_c, )
plot!(ylabel="Vial\nHeating\n", ylims=(0, 0.31), xlim=(0, NaN), legend=:right)
plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)

pl_mw = plot(xunit=u"hr", xwiden=false, ylabel="MW\nPower\n")
plot!(MWdat_c.t, MWdat_c.P, c=:red, label="", lw=0, fillrange=[0])
plot!(ylim=(0, 105))
plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)
pl_p = plot(xlabel="Time", xunit=u"hr", xwiden=false, ylabel=L"$p_\textrm{ch}$"*"\n")
plot!(pdat_c.t, pdat_c.pch, c=3, label="")
plot!(ylims=(0, 17), unitformat=:square)
plot!(xticks=(0:2:10), bottom_margin=20Plots.px)
pl_1c = plot(pl_T, pl_q, pl_mw, pl_p; size=(800,400), layout=@layout([a [b; c; d]]), link=:x, left_margin=20Plots.px)
end
savefig(plotsdir("gitter2019_1c.svg"))
savefig(plotsdir("gitter2019_1c.pdf"))

# ------------------ Figure 1c, ISLFD East Coast 2025 slides

begin
pl_T = blankplot_hrC(xwiden=false)
exptfplot!(T1dat_c.t, T1dat_c.T, T2dat_c.T, nmarks=40, )
plot!(Tsh, tmax=8u"hr", c=:black, label=L"T_\textrm{sh}")
modrftplot!(sol_c, labsuffix=", LC model")
# modrftplot!(sol_c, label=["f, LC-DIF" "vw, LC-DIF"])
# plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)

plot!(t3, T3, c=:purple, alpha=0.8, label=L"$T_\textrm{f}$, Park 2021")
plot!(legend=:topleft, )
plot!(bottom_margin=20Plots.px)

pl_q = qplotrf(sol_c, )
plot!(ylabel="Vial Heating\n", ylims=(0, 0.31), xlim=(0, NaN), legend=:right)
plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)

pl_mw = plot(xunit=u"hr", xwiden=false, ylabel="MW Power\n")
plot!(MWdat_c.t, MWdat_c.P, c=:red,lw=0, fillrange=[0], label="")
plot!(ylim=(0, 105))
plot!(xlabel="", xticks=(0:2:10, []), bottom_margin=-10Plots.px)
pl_p = plot(xlabel="Time", xunit=u"hr", xwiden=false, ylabel=L"Pirani $p_\textrm{ch}$"*"\n")
plot!(pdat_c.t, pdat_c.pch, c=3, label="")
plot!(ylims=(0, 17), unitformat=:square)
plot!(bottom_margin=20Plots.px)
pl_1c = plot(pl_T, pl_q, pl_mw, pl_p; size=(900,500), layout=@layout([a [b; c; d]]), link=:x)
plot!(left_margin=20Plots.px)
end
savefig(plotsdir("gitter2019_1c_slide.svg"))
savefig(plotsdir("gitter2019_1c_slide.pdf"))

# --------------------- All three in one big ugly figure

plot(pl_1a, pl_1b, pl_1c; layout=(1,3), size=(1200,600), left_margin=30Plots.px)
plot!(bottom_margin=20Plots.px)

plot(qplotrf(sol_a), qplotrf(sol_b), qplotrf(sol_c), layout=3, legend=(1.7, 0.7))
savefig(plotsdir("gitter2019_qplots.png"))
# ------------------
# Make fits into a table

# Get Bhambhani 2021 fit parameters
fit_bh_load = load(datadir("exp_pro", "bhambhani2021_fit_params.jld2"))
fit_bh = merge(fit_bh_load["conv"], fit_bh_load["mw"])

fits_gitter = transform.([trans_KKBBRp], [opt_a.u, opt_b.u, opt_c.u])
fit_table = Table(map([fit_bh, fits_gitter...]) do row
    (   Kshf = row.Kshf.val,
        Kvwf = row.Kvwf,
        Bf = row.Bf,
        Bvw = row.Bvw,
        a0 = row.Rp.R0,
        a1 = row.Rp.A1,
        a2 = row.Rp.A2)
end)


# For LaTeX output
# This \tabular gets copied and pasted into the .tex manuscript
# but inside the table with caption, \sisetup, etc.
# Also, replace {rllll} with {rSSSS} for siunitx formatting
# begin
# tab_transpose = hcat(collect.(fit_table)...)
# column_labels = LatexCell.(vcat(["\\makecell{Bhambhani 2021\\\\ Fig. 5b}"], ["\\makecell{Gitter 2019\\\\ Fig. 1"*i*"}" for i in "abc"]))
# row_labels = ["\$$(string(key)[1])\\s{$(string(key)[2:end])}\$, "*latexify(unit(val), fmt=SiunitxNumberFormatter()) for (key, val) in pairs(fit_table[1])]
# # table_formatter = (x,i,j)->latexify(ustrip(x), fmt=SiunitxNumberFormatter(format_options="round-mode=figures, round-precision=3, exponent-mode=threshold, exponent-thresholds=-2:3"))
# table_formatter = (x,i,j)->latexify(ustrip(x), fmt=SiunitxNumberFormatter())
# # table_highlighter = LatexHighlighter((x,i,j)->(i∈5:7 && j==1), (s,x,i,j)->s*"{\\footnotemark[1]}") # Footnote on Rp parameters
# pretty_table(tab_transpose; row_labels, column_labels, 
#     backend=:latex, 
#     alignment=:l,
#     formatters=[table_formatter],
#     # highlighters=[table_highlighter],
#     )
# end

# Same info, somewhat better LaTeX tabular output
# This gets put straight on the clipboard
# so you can paste it into the .tex manuscript
begin
column_labels = LaTeXString.(vcat(["{\\makecell{Bhambhani 2021\\\\ Fig. 5b}}"], ["{\\makecell{Gitter 2019\\\\ Fig. 1"*i*"}}" for i in "abc"]))
row_labels = LaTeXString.(["\$$(string(key)[1])\\s{$(string(key)[2:end])}\$, "*latexify(unit(val), fmt=SiunitxNumberFormatter()) for (key, val) in pairs(fit_table[1])])
clipboard(latextabular(ustrip.(tab_transpose),
    side=row_labels, 
    adjustment=[Symbol("@{}"), :r, :S, :S, :S, :S, Symbol("@{}")],
    head=column_labels,
    fmt=SiunitxNumberFormatter(),
    latex=false,
    booktabs=true))
end

# # -----------------
# # Braatz group analytical model

# const AM = MicrowaveLyoModeling.AnalyticalModel
# # They assume a sublimation temperature of 256.15K
# # This is equivalent to chamber pressure of 1300 μbar, since they do no mass transfer resistance
# # I will instead use measured pressures to get temperature
# Tm1 = find_zero( T->LyoPronto.calc_psub(T*u"K")-5u"μbar", 250)*u"K" # Lower bound on pch
# Tm2 = find_zero( T->LyoPronto.calc_psub(T*u"K")-20u"μbar", 250)*u"K" # Upper bound on pch
# Tm3 = 256.15u"K" # Their assumed value
# Kv = 65u"W/m^2/K"  # from paper
# Qppp_f = 242_345u"W/m^3"  # from paper
# kf = 2.30u"W/m/K"
# Tb0 = 236.85u"K" # initial T from Gitter 2019 temperatures
# Tbf = Tm3 # from Park 2021
# # Tb0 = 213.15u"K" # Gitter 2019 shelf temperatures
# # Tbf = 248.15u"K" # Gitter 2019 shelf temperatures for lyo
# r = 0.2u"K/minute"

# pp1 = AM.Params(Kv, Qppp_f, kf, 0.042u"m", r, Tb0, Tbf, Tm1)
# pp2 = AM.Params(Kv, Qppp_f, kf, 0.042u"m", r, Tb0, Tbf, Tm2)
# t1, s, T1 = AM.calc_tsT(pp1);
# t2, s, T2 = AM.calc_tsT(pp2);

# begin
# pl_all = blankplot_hrC()
# @df Tdat_b exptfplot!(:t, :T, nmarks=40)
# # @df thm_pd exptvwplot!(:t, :T3, nmarks=40)
# modrftplot!(sol_b, labsuffix=", LC-DIF",)
# plot!(t1, T1, c=:green, label=L"$T_{f}$, TLM with $T_{sub}$ from $p_{ch}$")
# plot!(t2, T2, c=:green, label="")
# plot!(t3, T3, c=:purple, label=L"$T_{f}$, TLM as originally shown")
# plot!(legend=:topleft)
# # plot!([1, 2.7], [-37, -27], arrow=:arrow, c=:gray, linewidth=2, label="")
# # pl_brtz = plot!(pl_all, u"hr", u"°C"; inset=bbox(0.38, 0.65, 0.3, 0.15), subplot=2)
# # pl_brtz = pl_all[2]
# # plot!(pl_brtz, t, T ,label="",  markersize=7, c=:green)
# # plot!(pl_brtz; ylim=(-40.1, -39.4), xlabel="", ylabel="", )
# end
# savefig(plotsdir("gitter2019_compT.svg"))
# savefig(plotsdir("gitter2019_compT.pdf"))
