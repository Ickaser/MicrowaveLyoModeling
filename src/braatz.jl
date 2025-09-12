"""
Implementation of the approximate analytical solution to the "hybrid freeze drying" model given in:

Prakitr Srisuma, George Barbastathis, Richard D. Braatz, "Analytical solutions for the modeling, 
optimization, and control of microwave-assisted freeze drying," Computers & Chemical Engineering,
Volume 177, 2023, 108318, ISSN 0098-1354, https://doi.org/10.1016/j.compchemeng.2023.108318.
"""
module AnalyticalModel

using Roots
using Unitful
using LyoPronto

struct Params
    Kv
    Qppp_f
    kf
    hf0
    r
    Tb0
    Tbf
    Tm
    ρCpf
    ρ
    ρa

    α
    σ
    λ
    ν
    κ1
    κ2
    Θ0
end

const ρ = 917u"kg/m^3"
const ρa = 63u"kg/m^3"
const ρCpf = ρ*1967.8u"J/kg/K"
const pice = 0.96

function Params(Kv, Qppp_f, kf, hf0, r, Tb0, Tbf, Tm)
    return Params(Kv, Qppp_f, kf, hf0, r, Tb0, Tbf, Tm, ρCpf, ρ, ρa, pice)
end

function Params(Kv, Qppp_f, kf, hf0, r, Tb0, Tbf, Tm, ρCpf, ρ, ρa, pice)
    α = kf/ρCpf
    σ = uconvert(NoUnits, r*hf0^2/(α*(Tm-Tb0)))
    λ = uconvert(NoUnits, Qppp_f*hf0/(Kv*(Tm-Tb0)))
    ν = uconvert(NoUnits, Kv*hf0/kf)
    κ1 = uconvert(NoUnits, r*Kv*hf0^3 / (α^2 * LyoPronto.ΔHsub*(ρ-ρa)*pice))
    κ2 = uconvert(NoUnits, hf0*(Kv*(Tb0-Tm) + Qppp_f*hf0)/(α*LyoPronto.ΔHsub*pice*(ρ-ρa)))
    Θ0 = 0.0

    return Params(Kv, Qppp_f, kf, hf0, r, Tb0, Tbf, Tm, ρCpf, ρ, ρa,
                    α, σ, λ, ν, κ1, κ2, Θ0)
                
end

τ(t, p::Params) = t/p.hf0^2*p.α
t(τ, p::Params) = τ*p.hf0^2/p.α

η(s, p::Params) = s/p.hf0
s(η, p::Params) = η*p.hf0

Θ(T, p::Params) = (T-p.Tb0)/(p.Tm-p.Tb0)
T(Θ, p::Params) = Θ*(p.Tm-p.Tb0) + p.Tb0
# Hw is volumetric heat generation in water
# pbw is fraction water

# Hv = Hw*pbw volumetric heat generation in frozen domain
# Tb = Tb0 + rt # shelf ramp
# Hb = h*(T(L) - Tb)

# s: 0 at start, goes to L
# x=0 top, x=L bottom

# η = s/L
# τ = αt/L^2, α = k/ρCp
# ν = Kv*L/k # Biot number in ice
# σ = r*L^2/α/(Tm-Tb0)
#Θ = (T - Tb0)/(Tm - Tb0)

# ξ = x/L
# τm: θ = θm
# τmax: 

function calc_Θ_heating(τ, ξ, p::Params)
    term1 = ((p.λ/2 - p.σ/2 - 3*p.ν*p.Θ0/(2p.ν + 6))*exp(-3p.ν/(p.ν+3)*τ) - p.λ/2 + p.σ/2)
    Θ = p.σ*τ-(1 + 2/p.ν)*term1 + ξ^2*term1
    return Θ
end

function calc_tsT(p::Params)

    τm = find_zero(t->calc_Θ_heating(t, 0, p), 1)
    τmax = τ((p.Tbf - p.Tb0)/p.r, p)

    τ1 = range(0.0, τm, length=15)
    Θ1 = calc_Θ_heating.(τ1, 1, (p,))
    η1 = fill(0.0, length(τ1))

    τ2 = range(τm, τmax, length=15)[begin+1:end]
    η2 = @. p.κ1*(τ2^2-τm^2)/2 + p.κ2*(τ2-τm)

    ηmax = η2[end]
    τf = (1.0 - ηmax)/(p.κ1*τmax + p.κ2) + τmax
    τ3 = range(τmax, τf, length=15)[begin+1:end]
    η3 = range(ηmax, 1.0, length=15)[begin+1:end]

    Θ2 = fill(1.0, length(τ2))
    Θ3 = fill(1.0, length(τ3))

    τr = vcat(τ1, τ2, τ3)
    tr = uconvert.(u"hr", t.(τr, (p,)))
    Θr = vcat(Θ1, Θ2, Θ3)
    Tr = uconvert.(u"K", T.(Θr, (p,)))
    ηr = vcat(η1, η2, η3)
    sr = uconvert.(u"cm", s.(ηr, (p,)))
    return tr, sr, Tr
end


end # module AnalyticalModel