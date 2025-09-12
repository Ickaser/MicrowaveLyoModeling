export plot_defaults_lprf
function plot_defaults_lprf()
    default(:fontfamily, "Computer Modern")
    default(:framestyle, :box)
    default(:lw, 2)
    default(:markersize, 4)
    default(:markerstrokewidth, 0.5)
    resetfontsizes()
    scalefontsizes(1.2)
end

export qplotrf
function qplotrf(sol; tot_lab="Total", kw...)
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
    names = ["RF-f", "vw-f", "sh-f"]
    labs = ["\$Q_\\textrm{$nm}\$" for nm in names]
    pl = plot(u"hr", u"W", xlabel="Time", ylabel="Heating", unitformat=:square; kw...)
    areastackplot!(sol.t, QRFf, Qvwf, Qshf, labels=permutedims(labs), fillalpha=0.6)
    plot!(sol.t, Qshf.+Qvwf.+QRFf, c=:black, label=tot_lab)
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

