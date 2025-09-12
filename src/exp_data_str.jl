export lyostar_columnrename, revo_columnrename
export fiberoptic_columnrename, powermeter_columnrename
export save_fitresults
export identify_pd_end_2der, identify_pd_end_onoff
"""
Currently missing a way to handle any potentially-present thermocouple values.
"""
function lyostar_columnrename(row)
    # nt = (tstamp = DateTime(row.Timestamp, dateformat"mm/dd/yyyy H:M:S"),
    nt = (tstamp = row.Timestamp,
          pch_sp = row.var"SPLYO.VACUUM_SP.F_CV"*u"mTorr",
          pch_pir = row.var"SPLYO.CHAMBER_PIRANI.F_CV"*u"mTorr",
          pch_cm = row.var"SPLYO.CHAMBER_CM.F_CV"*u"mTorr",
          Tsh_sp = row.var"SPLYO.SHELF_SP.F_CV"*u"°C",
          Tsh_i = row.var"SPLYO.SHELF_INLET.F_CV"*u"°C",
          Tsh_o = row.var"SPLYO.SHELF_OUTLET.F_CV"*u"°C",
    )
    return nt
end

"""
Currently missing a way to handle any potentially-present thermocouple values.
"""
function revo_columnrename(row)
    # nt = (tstamp = DateTime(row.Timestamp, dateformat"mm/dd/yyyy H:M:S"),
    nt = (tstamp = row.CycleTime,
          phase = row.var"Phase",
          step = row.var"Step",
          pch_sp = row.var"V-SetPT"*u"mTorr",
          pch_pir = row.var"Vacuum-P"*u"mTorr",
          pch_cm = row.var"Vacuum-c"*u"mTorr",
          Tsh_sp = row.var"S-SetPT"*u"°C",
          Tsh_i = row.var"Shelf"*u"°C",
    )
    return nt
end


function fiberoptic_columnrename(row)
    nt = (tstamp = row.var"Elapsed [s]"*u"s",
          T1 = row.var"T1 [C]"*u"°C",
          T2 = row.var"T2 [C]"*u"°C",
          T3 = row.var"T3 [C]"*u"°C",
          T4 = row.var"T4 [C]"*u"°C",
    )
    return nt
end
function powermeter_columnrename(row)
    nt = (tstamp = row.var"Elapsed [s]"*u"s",
          f = row.var"Frequency [GHz]"*u"GHz",
          P = row.var"Power [W]"*u"W",
    )
    return nt
end

function save_fitresults(opt::Union{SciMLBase.OptimizationSolution, SciMLBase.NonlinearSolution}, casename; K_edge = nothing)
    pass = opt.prob.p
    # trans, po, fitdat =
    fitprm = transform(pass[1], opt.u)
    po = setproperties(pass[2], fitprm)
    save_fitresults(po, pass[3], casename; K_edge)
end

function save_fitresults(po::ParamObjRF, fitdat::PrimaryDryFit, casename; K_edge = nothing)
    sol = solve(ODEProblem(po), Rodas3())
    Tferr = obj_expT(sol, fitdat, tweight=0, Tvw_weight=0)
    Tvwerr = obj_expT(sol, fitdat, tweight=0, Tvw_weight=1) - Tferr
    terr = obj_expT(sol, fitdat, tweight=1, Tvw_weight=0) - Tferr
    Kshf = po.Kshf(po.pch(0u"s"))
    α = po.alpha
    (; Rp, Bf, Bvw, Kvwf) = po
    allfit = @strdict Kshf Rp K_edge Kvwf α  Bf Bvw po terr Tferr Tvwerr
    safesave(datadir("exp_pro", "$(casename)_KvRpRF.jld2"), allfit)
    heating = qrf_integrate(sol, po)
    safesave(datadir("exp_pro", "$(casename)_Q.jld2"), heating) 
end

function identify_pd_end_2der(data, window_width=91, tmin=0u"hr", tmax=Inf*u"hr")
    pch_pir_der2 = savitzky_golay(data.pch_pir, window_width, 3, deriv=2).y
    t_end = data.t[argmax(pch_pir_der2[tmin .< data.t .< tmax])] + tmin
    return t_end
end
"""
Identify a range of times for the end of primary drying, with the onset and offset of the
Pirani pressure curve.

Use a Savitzky-Golay filter to smooth, identify minimum of derivative, and then find linear
intersects with maximum and minimum values of pressure. Those intersects are taken as onset 
and offset.
"""
function identify_pd_end_onoff(data, window_width=91, tmin=0u"hr", tmax=Inf*u"hr")
    pch_pir_sm = savitzky_golay(data.pch_pir, window_width, 3, deriv=0).y
    dt = sum(diff(data.t)) / (length(data.t)-1)
    pch_pir_der1 = savitzky_golay(data.pch_pir, window_width, 3, deriv=1, rate=1/dt).y
    ti_inds = tmin .< data.t .< tmax
    i_mid = argmin(pch_pir_der1[ti_inds]) + searchsortedfirst(data.t, tmin)
    t_mid = data.t[i_mid]
    dp_mid = pch_pir_der1[i_mid]
    p_mid = pch_pir_sm[i_mid]
    p_max = maximum(pch_pir_sm[ti_inds])
    p_min = minimum(pch_pir_sm[ti_inds])
    t_onset = (p_max - p_mid)/dp_mid + t_mid
    t_offset = (p_min - p_mid)/dp_mid + t_mid
    return (t_onset, t_offset)
end