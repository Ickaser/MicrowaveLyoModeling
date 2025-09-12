export lumped_cap_rf_alt2, lumped_cap_rf_alt1
export gen_sol_rf_LC1, gen_sol_rf_LC2, gen_sol_rf_LC3, gen_sol_rf_LC3b
export obj_KBBa1, obj_KBBa2, obj_KBBaa3

const k_sucrose = LyoPronto.k_sucrose
const rho_glass = LyoPronto.rho_glass
const σ = LyoPronto.σ
const ΔHsub = LyoPronto.ΔHsub
const S_interp = LyoPronto.S_interp
const e_0 = LyoPronto.e_0

const rfprm_base_scale = (1.0u"W/m^2/K", 1e7u"Ω/m^2", 1e7u"Ω/m^2", 0.1u"cm^1.5")

function gen_sol_rf_LC3(KBB_log, po::ParamObjRF; rfprm_scale = rfprm_base_scale, kwargs...)
    newp = copy(po)
    @reset newp.K_vwf = exp(KBB_log[1])*rfprm_scale[1]
    @reset newp.B_f = exp(KBB_log[2])*rfprm_scale[2]
    @reset newp.B_vw = exp(KBB_log[3])*rfprm_scale[3]
    newprob = ODEProblem(newp)
    sol = solve(newprob, Rodas3(); kwargs...)
    return sol
end
function gen_sol_rf_LC3(KBB_log, po::ParamObjRF, fitdat::PrimaryDryFit; rfprm_scale = rfprm_base_scale, kwargs...)
    newp = copy(po)
    @reset newp.K_vwf = exp(KBB_log[1])*rfprm_scale[1]
    @reset newp.B_f = exp(KBB_log[2])*rfprm_scale[2]
    @reset newp.B_vw = exp(KBB_log[3])*rfprm_scale[3]
    newprob = ODEProblem(newp)
    sol = solve(newprob, Rodas3(); saveat=ustrip.(u"hr", fitdat.t), kwargs...)
    return sol
end

function gen_sol_rf_LC3b(KBBaa_log, po::ParamObjRF; rfprm_scale = rfprm_base_scale, kwargs...)
    newp = copy(po)
    @reset newp.K_vwf = exp(KBBaa_log[1])*rfprm_scale[1]
    @reset newp.B_f = exp(KBBaa_log[2])*rfprm_scale[2]
    @reset newp.B_vw = exp(KBBaa_log[3])*rfprm_scale[3]
    @reset newp.Rp = RpFormFit(newp.Rp.R0, exp(KBBaa_log[4])*10u"cm*Torr*hr/g", exp(KBBaa_log[5])*0.1u"cm^-1")
    newprob = ODEProblem(newp)
    sol = solve(newprob, Rodas3(); kwargs...)
    return sol
end
function gen_sol_rf_LC3b(KBBaa_log, po::ParamObjRF, fitdat::PrimaryDryFit; rfprm_scale = rfprm_base_scale, kwargs...)
    newp = copy(po)
    @reset newp.K_vwf = exp(KBBaa_log[1])*rfprm_scale[1]
    @reset newp.B_f = exp(KBBaa_log[2])*rfprm_scale[2]
    @reset newp.B_vw = exp(KBBaa_log[3])*rfprm_scale[3]
    @reset newp.Rp = RpFormFit(newp.Rp.R0, exp(KBBaa_log[4])*10u"cm*Torr*hr/g", exp(KBBaa_log[5])*0.1u"cm^-1")
    newprob = ODEProblem(newp)
    sol = solve(newprob, Rodas3(); saveat=ustrip.(u"hr", fitdat.t), kwargs...)
    return sol
end

# -----------------
# Old model forms
# Note that the sqrt(hd) term has a derivative of Inf at hd=0,
# so these models don't play nice with autodiff, 
# unless we were to perturb them away from 0 (which we could do...)

function gen_sol_rf_LC2(fitlog, tr, po::ParamObjRF; badprms=nothing, saveat=[], kwargs...)
    fitprm = transform(tr, fitlog)
    prms = setproperties(po, fitprm)
    !isnothing(badprms) && badprms(prms) && return NaN
    prob = ODEProblem(lumped_cap_rf_alt2, LyoPronto.calc_u0(prms), (0.0, 1000.0), prms)
    sol = solve(prob, Rosenbrock23(autodiff=AutoFiniteDiff()); saveat, callback=end_drying_callback, kwargs...)
end
function gen_sol_rf_LC2(fitlog, tr, po, fitdat; badprms=nothing, kwargs...)
    sol = gen_sol_rf_LC2(fitlog, tr, po; saveat=ustrip.(u"hr", fitdat.t), badprms, kwargs...)
    return sol
end

function gen_sol_rf_LC1(fitlog, tr, po; saveat=[], badprms=nothing, kwargs...)
    fitprm = transform(tr, fitlog)
    prms = setproperties(po, fitprm)
    !isnothing(badprms) && badprms(prms) && return NaN
    prob = ODEProblem(lumped_cap_rf_alt1, LyoPronto.calc_u0(prms), (0.0, 1000.0), prms)
    sol = solve(prob, Rosenbrock23(autodiff=AutoFiniteDiff()); saveat, callback=end_drying_callback, kwargs...)
    return sol
end
function gen_sol_rf_LC1(fitlog, tr, po, fitdat; badprms=nothing, kwargs...)
    sol = gen_sol_rf_LC1(fitlog, tr, po; saveat=ustrip.(u"hr", fitdat.t), badprms, kwargs...)
    return sol
end

@doc raw"""
    lumped_cap_rf_alt2(u, params, tn)

Compute the right-hand-side function for the ODEs making up the lumped-capacitance microwave-assisted model.

An alternate model to [`lumped_cap_rf`](@ref); syntax shared there.
"""
function lumped_cap_rf_alt2(du, u, params, tn)
    du .= lumped_cap_rf_LC2(u, params, tn)[1]
end

@doc raw"""
    lumped_cap_rf_LC2(u, params, tn)

This does the work for [`lumped_cap_rf_alt2`](@ref), but returns `dudt,  [Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw]` with `Q_...` as Unitful quantities in watts. 
The extra results are helpful in investigating the significance of the various heat transfer modes,
but are not necessary in the ODE integration.

LC2: Q_shw evaluated with radiation; shape factor included; α>=0
The version used for NIIMBL final report.
"""
function lumped_cap_rf_LC2(u, params, tn)
    # Unpack all the parameters
    Rp, h_f0, cSolid, ρ_solution = params[1]
    K_shf_f, A_v, A_p, = params[2]
    pch, Tsh, P_per_vial = params[3] 
    m_f0, cp_f, m_v, cp_v, A_rad = params[4]
    f_RF, epp_f, epp_vw = params[5]
    K_vwf, B_f, B_vw, alpha = params[6]
    # Dimensionalize the state variables
    t = tn*u"hr" 
    m_f = u[1]*u"g"
    T_f = u[2]*u"K"
    T_vw = u[3]*u"K"
    # Compute some properties
    porosity = (ρ_solution - cSolid)/ρ_solution
    k_dry = k_sucrose*(1-porosity)
    V_vial = m_v / rho_glass
    # Do some geometry
    rad = sqrt(A_p/π)
    h_f = m_f/m_f0 * h_f0 
    h_d = h_f0 - h_f
    # Heat transfer from shelf
    K_shf = K_shf_f(pch(t))
    Q_shf = K_shf*A_v*(Tsh(t)-T_f) 
    Q_shw = 0.9*σ*(Tsh(t)^4-T_vw^4)*A_rad # As if exchanging heat above vial
    # Evaluate mass flow; positive means drying is progressing. Not forced to be positive
    A_sub = A_p + alpha*sqrt(abs(h_d+ .00001u"cm"))
    mflow = A_sub/Rp(h_d)*(calc_psub(T_f) - pch(t)) # g/s
    Q_sub = mflow*ΔHsub # Sublimation
    # Evaluate heat transfer from wall
    Bi = uconvert(NoUnits, K_vwf*rad/k_dry)
    Q_vwf = 2π*(K_vwf*rad*h_f + k_dry*(h_f0-h_f)*S_interp(Bi)) * (T_vw-T_f)
    # Volumetric heating
    Qppp_RF_f  = 2*pi*f_RF*e_0*epp_f(T_f, f_RF)*P_per_vial(t)*B_f # W / m^3
    Qppp_RF_vw = 2*pi*f_RF*e_0*epp_vw*P_per_vial(t)*B_vw # W / m^3
    Q_RF_f = Qppp_RF_f*A_p*h_f # W
    Q_RF_vw = Qppp_RF_vw*V_vial # W
    # Check that total volumetric heating is less than input power
    if Q_RF_f + Q_RF_vw > P_per_vial(t) && t == 0u"hr"
        @warn "Energy balance of EM terms not satisfied." uconvert(u"W", Q_RF_f) uconvert(u"W", Q_RF_vw) P_per_vial(t)
    end
    # Evaluate derivatives
    # Desublimation is not allowed here: if we clamp mflow itself, then the DAE is unstable
    dm_f = min(0.0u"kg/s", -mflow/porosity)
    dT_f = (Q_RF_f-Q_sub+Q_vwf+Q_shf) / (m_f*cp_f) - T_f*dm_f/m_f
    dT_vw = (Q_RF_vw+Q_shw-Q_vwf) / (m_v*cp_v)
    # Strip units from derivatives; return all heat transfer terms
    return ustrip.([u"g/hr", u"K/hr", u"K/hr"], [dm_f, dT_f, dT_vw]), uconvert.(u"W", [Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw, ])
end

@doc raw"""
    lumped_cap_rf_alt1(u, params, tn)

Compute the right-hand-side function for the ODEs making up the lumped-capacitance microwave-assisted model.

An alternate model to [`lumped_cap_rf`](@ref); syntax shared there.
"""
function lumped_cap_rf_alt1(du, u, params, tn)
    du .= lumped_cap_rf_LC1(u, params, tn)[1]
end

@doc raw"""
    lumped_cap_rf_LC1(u, params, tn)

This does the work for [`lumped_cap_rf_alt3`](@ref), but returns `dudt,  [Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw]` with `Q_...` as Unitful quantities in watts. 
The extra results are helpful in investigating the significance of the various heat transfer modes,
but are not necessary in the ODE integration.
 
LC1: Q_shw evaluated with radiation; no shape factor; α>=0
The version used for ScientificReports article.
"""
function lumped_cap_rf_LC1(u, params, tn)
    # Unpack all the parameters
    Rp, h_f0, cSolid, ρ_solution = params[1]
    K_shf_f, A_v, A_p, = params[2]
    pch, Tsh, P_per_vial = params[3] 
    m_f0, cp_f, m_v, cp_v, A_rad = params[4]
    f_RF, epp_f, epp_vw = params[5]
    K_vwf, B_f, B_vw, alpha = params[6]
    # Dimensionalize the state variables
    t = tn*u"hr" 
    m_f = u[1]*u"g"
    T_f = u[2]*u"K"
    T_vw = u[3]*u"K"
    # Compute some properties
    porosity = (ρ_solution - cSolid)/ρ_solution
    k_dry = k_sucrose*(1-porosity)
    V_vial = m_v / rho_glass
    # Do some geometry
    rad = sqrt(A_p/π)
    h_f = m_f/m_f0 * h_f0 
    h_d = h_f0 - h_f
    # Heat transfer from shelf
    K_shf = K_shf_f(pch(t))
    Q_shf = K_shf*A_v*(Tsh(t)-T_f) 
    Q_shw = 0.9*σ*(Tsh(t)^4-T_vw^4)*A_rad # As if exchanging heat above vial
    # Evaluate mass flow; positive means drying is progressing. Not forced to be positive
    A_sub = A_p + alpha*sqrt(abs(h_d + .00001u"cm"))
    mflow = A_sub/Rp(h_d)*(calc_psub(T_f) - pch(t)) # g/s
    Q_sub = mflow*ΔHsub # Sublimation
    # Evaluate heat transfer from wall
    Q_vwf = 2π*K_vwf*rad*h_f*(T_vw-T_f)
    # Volumetric heating
    Qppp_RF_f  = 2*pi*f_RF*e_0*epp_f(T_f, f_RF)*P_per_vial(t)*B_f # W / m^3
    Qppp_RF_vw = 2*pi*f_RF*e_0*epp_vw*P_per_vial(t)*B_vw # W / m^3
    Q_RF_f = Qppp_RF_f*A_p*h_f # W
    Q_RF_vw = Qppp_RF_vw*V_vial # W
    # Check that total volumetric heating is less than input power
    if Q_RF_f + Q_RF_vw > P_per_vial(t) && t == 0u"hr"
        @warn "Energy balance of EM terms not satisfied." uconvert(u"W", Q_RF_f) uconvert(u"W", Q_RF_vw) P_per_vial(t)
    end
    # Evaluate derivatives
    # Desublimation is not allowed here: if we clamp mflow itself, then the DAE is unstable
    dm_f = min(0u"kg/s", -mflow/porosity)
    dT_f = (Q_RF_f-Q_sub+Q_vwf+Q_shf) / (m_f*cp_f) - T_f*dm_f/m_f
    dT_vw = (Q_RF_vw+Q_shw-Q_vwf) / (m_v*cp_v)
    # Strip units from derivatives; return all heat transfer terms
    return ustrip.([u"g/hr", u"K/hr", u"K/hr"], [dm_f, dT_f, dT_vw]), uconvert.(u"W", [Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw, ])
end
