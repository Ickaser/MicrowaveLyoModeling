export plot_defaults_mlm
function plot_defaults_mlm()
    default(:fontfamily, "Computer Modern")
    default(:framestyle, :box)
    default(:lw, 2)
    default(:markersize, 4)
    default(:markerstrokewidth, 0.5)
    default(:unitformat, :square)
    resetfontsizes()
    scalefontsizes(1.2)
    set_default(labelformat=:square) # Latexify, not Plots
end

export qplotrf
function qplotrf(sol; tot_lab="Total", ordering=1:3, kw...)
    RF_params = sol.prob.p
    Qcontrib = map(sol.t) do ti
        LyoPronto.lumped_cap_rf!(fill(0.0, 3), sol(ti), RF_params, ti, Val(true))
    end
    Qcontrib = hcat(Qcontrib...)
    Qsub = Qcontrib[1,:]
    Qshf = Qcontrib[2,:]
    Qvwf = Qcontrib[3,:]
    QRFf = Qcontrib[4,:]
    QRFvw = Qcontrib[5,:]
    # names = ["sub", "sh-f", "vw-f", "RF-f", "RF-vw", "sh-vw"]
    modes = [QRFf, Qvwf, Qshf][ordering]
    names = ["RF-f", "vw-f", "sh-f"][ordering]
    labs = ["\$Q_\\textrm{$nm}\$" for nm in names]
    pl = plot(u"hr", u"W", xlabel="Time", ylabel="Heating", unitformat=:square; )
    fillcolor = permutedims([:yellow, :orange, :red][ordering])
    plot!(sol.t, Qshf.+Qvwf.+QRFf, c=:black, label=tot_lab)
    areastackplot!(sol.t, modes...; labels=permutedims(labs), fillalpha=0.6, fillcolor, kw...)
    plot!(xlabel="Time", ylabel="Heating")
    return pl
end
export qplotrf_alt2
function qplotrf_alt2(sol; kw...)
    RF_params = sol.prob.p
    Qcontrib = map(sol.t) do ti
        lumped_cap_rf_LC2(sol(ti), RF_params, ti)[2]
    end
    Qcontrib = hcat(Qcontrib...)
    Qsub = Qcontrib[1,:]
    Qshf = Qcontrib[2,:]
    Qvwf = Qcontrib[3,:]
    QRFf = Qcontrib[4,:]
    QRFvw = Qcontrib[5,:]
    # names = ["sub", "sh-f", "vw-f", "RF-f", "RF-vw", "sh-vw"]
    names = ["RF-f", "vw-f", "sh-f"]
    labs = ["\$Q_\\text{$nm}\$" for nm in names]
    pl = plot(u"hr", u"W", xlabel="Time", ylabel="Heating", unitformat=:square; kw...)
    areastackplot!(sol.t, QRFf, Qvwf, Qshf, labels=permutedims(labs))
    plot!(sol.t, Qshf.+Qvwf.+QRFf, c=:black, label="Total")
    return pl
end

export qplotrf_alt1
function qplotrf_alt1(sol; kw...)
    RF_params = sol.prob.p
    Qcontrib = map(sol.t) do ti
        lumped_cap_rf_LC1(sol(ti), RF_params, ti)[2]
    end
    Qcontrib = hcat(Qcontrib...)
    Qsub = Qcontrib[1,:]
    Qshf = Qcontrib[2,:]
    Qvwf = Qcontrib[3,:]
    QRFf = Qcontrib[4,:]
    QRFvw = Qcontrib[5,:]
    # names = ["sub", "sh-f", "vw-f", "RF-f", "RF-vw", "sh-vw"]
    names = ["RF-f", "vw-f", "sh-f"]
    labs = ["\$Q_\\text{$nm}\$" for nm in names]
    pl = plot(u"hr", u"W", xlabel="Time", ylabel="Heating"; kw...)
    areastackplot!(sol.t, QRFf, Qvwf, Qshf, labels=permutedims(labs))
    plot!(sol.t, Qshf.+Qvwf.+QRFf, c=:black, label="Total")
    return pl
end



export blankplot_hrC
function blankplot_hrC(;kwargs...)
    plot(u"hr", u"°C", xlabel="Time", ylabel="Temperature", unitformat=:square; kwargs...)
end

