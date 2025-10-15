# -------------
# Notes:
# This file provides postprocessing of the experimental data used in this project.
# The below functions are only intended to be run once (assuming nothing needs 
# to be changed); they are encapsulated in `let` blocks to 
# isolate each experiment into its own scope.

# ----------------------------
# M1 

let 
dat1 = CSV.read(datadir("exp_raw", "M1", "TemperatureLog_030223.csv"), Table)
dat2 = CSV.read(datadir("exp_raw", "M1", "Mannitol_5ml_Temp.csv"), Table, 
    comment="#", stripwhitespace=true, dateformat="mm/dd/yyyy H:M:S")

time1 = dat1.var"Elapsed [s]"

lyo_full_pre = map(lyostar_columnrename, dat2)
time2 = uconvert.(u"hr", lyo_full_pre.tstamp .- lyo_full_pre.tstamp[1])
lyo_full = Table(lyo_full_pre; t=time2)

# Conditions for identifying start and end of primary drying data
istart_lyo_pd = findfirst(lyo_full.pch_sp .== 100u"mTorr")
iend_lyo_pd = length(lyo_full)

lyo_pd = Table(lyo_full[istart_lyo_pd:iend_lyo_pd])
lyo_pd.t .-= lyo_pd.t[1]

# Plot temperature and pressure to make sure things look right
# plT = plot(lyo_pd.Tsh_sp)
# plp = plot(lyo_pd.pch_pir, ylim=(0, 300))
# plot!(plp, lyo_pd.pch_cm)
# plot(plT, plp, layout=(2,1))

# Continue reading data
thm_full = map(fiberoptic_columnrename, dat1)

# Conditions for identifying start and end of primary drying
istart_thm_pd = argmin(thm_full.T1)
iend_thm_pd = argmax(thm_full.T4) + 600

thm_pd_sub = Table(thm_full[istart_thm_pd:iend_thm_pd])
thm_pd = Table(thm_pd_sub; t = uconvert.(u"hr", thm_pd_sub.tstamp .- thm_pd_sub.tstamp[1]))

# Identify second derivative condition
pch_pir_sm = savitzky_golay(lyo_pd.pch_pir, 91, 3, deriv=0).y 
pch_pir_d2 = savitzky_golay(lyo_pd.pch_pir, 91, 3, rate = 1u"minute", deriv=2).y
t_end = lyo_pd.t[argmax(pch_pir_d2[100:end])+99] - lyo_pd.t[1] +0.0u"hr" # Make it be a float

# Sanity check for second derivative
# @df lyo_pd plot(:t, :pch_pir)
# @df plot!(:t, pch_pir_sm)
# tendplot!(t_end)

Tsh_rv = RampedVariable(uconvert.(u"K", [-40u"°C", 10u"°C"]).+0.0u"K", [0.5u"K/minute"], [])

# Plot temperatures for help in identifying type and useful region
# labels = ["v1, T4" "v2, T1" "gl, T3"]
# @df thm_pd exptfplot(:t, :T4, :T1, :T3, label=labels)
# plot!(Tsh_rv, label="shelf")

iend_T4 = argmin(thm_pd.T4[1000:end]) + 1000
iend_T1 = argmin(thm_pd.T1[1000:end]) + 1000

# Condense down into the part that will be fit
fitdat = @df thm_pd[begin:6:end] PrimaryDryFit(:t, (:T4[begin:iend_T4÷6], 
    :T1[begin:iend_T1÷6]), :T3, t_end);

# Save data in its postprocessed form, to avoid repeating the above
safesave(datadir("exp_pro", "M1_processed.jld2"), @strdict lyo_pd thm_pd thm_full lyo_full t_end fitdat)
end

# ----------------------------
# M2 

let 
dat1 = CSV.read(datadir("exp_raw", "M2", "TemperatureLog_021023.csv"), Table)
dat2 = CSV.read(datadir("exp_raw", "M2", "RF_Petr_5_M_TEMP.csv"), Table, 
    comment=";", stripwhitespace=true, dateformat="mm/dd/yyyy H:M:S")

time1 = dat1.var"Elapsed [s]"

lyo_full_pre = map(lyostar_columnrename, dat2)
time2 = uconvert.(u"hr", lyo_full_pre.tstamp .- lyo_full_pre.tstamp[1])
lyo_full = Table(lyo_full_pre; t=time2)

# Conditions for identifying start and end of primary drying data
istart_lyo_pd = findfirst(lyo_full.pch_sp .== 100u"mTorr")
iend_lyo_pd = length(lyo_full)

lyo_pd = Table(lyo_full[istart_lyo_pd:iend_lyo_pd])
lyo_pd.t .-= lyo_pd.t[1]

# Plot temperature and pressure to make sure things look right
# plT = plot(lyo_pd.Tsh_sp)
# plp = plot(lyo_pd.pch_pir, ylim=(0, 300))
# plot!(plp, lyo_pd.pch_cm)
# plot(plT, plp, layout=(2,1))

# Continue reading data
thm_full_pre = map(fiberoptic_columnrename, dat1)
thm_full = Table(thm_full_pre; t = uconvert.(u"hr", thm_full_pre.tstamp))

# Conditions for identifying start and end of primary drying
istart_thm_pd = argmin(thm_full.T1)+204
iend_thm_pd = length(thm_full)

thm_pd_sub = Table(thm_full[istart_thm_pd:iend_thm_pd])
thm_pd = Table(thm_pd_sub; t = uconvert.(u"hr", thm_pd_sub.tstamp .- thm_pd_sub.tstamp[1]))

# Identify second derivative condition
pch_pir_sm = savitzky_golay(lyo_pd.pch_pir, 91, 3, deriv=0).y
pch_pir_d2 = savitzky_golay(lyo_pd.pch_pir, 91, 3, deriv=2, rate=1u"minute").y
t_end = lyo_pd.t[argmax(pch_pir_d2[100:end])+99] - lyo_pd.t[1] +0.0u"hr"

# Sanity check for second derivative
# @df lyo_pd plot(:t, :pch_pir)
# @df plot!(:t, pch_pir_sm)
# tendplot!(t_end)

Tsh_rv = RampedVariable(uconvert.(u"K", [-40u"°C", 10u"°C"]).+0.0u"K", [0.5u"K/minute"], [])

# Plot temperatures for help in identifying type and useful region
@df thm_pd exptfplot(:t, :T1, :T3, :T4)
plot!(Tsh_rv, label="shelf", tmax=10u"hr")

iend_T1 = searchsortedfirst(thm_pd.t, 3u"hr")
iend_T34 = searchsortedfirst(thm_pd.t, 5u"hr")

# Condense down into the part that will be fit
fitdat = @df thm_pd PrimaryDryFit(:t, (:T4[begin:iend_T34], 
    :T1[begin:iend_T1]), maximum(:T4), t_end);

# Save data in its postprocessed form, to avoid repeating the above
safesave(datadir("exp_pro", "M2_processed.jld2"), @strdict lyo_pd thm_pd thm_full lyo_full t_end fitdat)
end


# ------------------------
# M3

let 
dat1 = CSV.read(datadir("exp_raw", "M3", "CL_2_21_5M.csv"), Table, 
    dateformat="mm/dd/yyyy H:M",comment=";", stripwhitespace=true)
dat2 = CSV.read(datadir("exp_raw", "M3", "FiberOptic_21Feb2023.csv"), Table, comment="#")
dat3 = CSV.read(datadir("exp_raw", "M3", "PowerMeter_21Feb2023.csv"), Table, comment="#", stripwhitespace=true)

time2 = dat2.var"Elapsed [s]"

lyo_full_pre = map(lyostar_columnrename, dat1)
time1 = uconvert.(u"hr", lyo_full_pre.tstamp .- lyo_full_pre.tstamp[1])
lyo_full = Table(lyo_full_pre; t=time1)

# Conditions for identifying start and end of primary drying data
istart_lyo_pd = findfirst(lyo_full.Tsh_sp[250:end] .== -15u"°C") + 250
iend_lyo_pd = findlast(lyo_full.pch_pir .> 95u"mTorr")

lyo_pd = Table(lyo_full[istart_lyo_pd:iend_lyo_pd])
lyo_pd.t .-= lyo_pd.t[1]

# Plot temperature and pressure to make sure things look right
plT = plot(lyo_pd.Tsh_sp)
plp = plot(lyo_pd.pch_pir, ylim=(0, 300))
plot!(plp, lyo_pd.pch_cm)
plot(plT, plp, layout=(2,1))

# Identify second derivative condition
pch_pir_sm = savitzky_golay(lyo_pd.pch_pir, 91, 3, deriv=0).y 
pch_pir_d2 = savitzky_golay(lyo_pd.pch_pir, 91, 3, rate=1u"minute", deriv=2).y
t_end = lyo_pd.t[argmax(pch_pir_d2[100:end])+99] - lyo_pd.t[1] +0.0u"hr"

# Sanity check for second derivative
# @df lyo_pd plot(:t, :pch_pir)
# @df plot!(:t, pch_pir_sm)
# tendplot!(t_end)

# Continue reading data
thm_full = map(fiberoptic_columnrename, dat2)

# Conditions for identifying start and end of primary drying
istart_thm_pd = findlast(thm_full.tstamp .< 10u"s")
# iend_thm_pd = findlast(thm_full.tstamp.-thm_full.tstamp[istart_thm_pd] .<10u"hr")
iend_thm_pd = length(thm_full)

thm_pd_sub = Table(thm_full[istart_thm_pd:iend_thm_pd])
thm_pd = Table(thm_pd_sub; t = uconvert.(u"hr", thm_pd_sub.tstamp .- thm_pd_sub.tstamp[1]))

# Plot temperatures for help in identifying type and useful region
labels = ["v1, T4" "v2, T3" "gl, T1"]
@df thm_pd exptfplot(:t, :T4, :T3)
@df thm_pd exptvwplot!(:t, :T1, trim=10)
@df lyo_pd plot!(:t, :Tsh_sp)


# Condense down into the part that will be fit
vw_valid = thm_pd.T1 .> -10u"°C"
fitdat = @df thm_pd[begin:12:end] PrimaryDryFit(:t, (:T4, 
    :T3[:T3 .< -5u"°C"]), :T1, t_end);

# Special: power output from power meter
pow_full = map(powermeter_columnrename, dat3)
istart_pow = findlast(pow_full.tstamp .<= 5u"s")
pow_pd = pow_full[istart_pow:end]
P_per_vial = LinearInterpolation(pow_pd.P .*0.54./17, pow_pd.tstamp, extrapolation=ExtrapolationType.Linear) # Correct for efficiency, # of vials

safesave(datadir("exp_pro", "M3_processed.jld2"), @strdict lyo_pd thm_pd thm_full lyo_full t_end fitdat P_per_vial)

end

# -------------------------
# BSA conventional

let
dat1 = CSV.read(datadir("exp_raw", "BSA", "BSA_Single_Shelf_CONV", "2023-09-12-13_REVO_AD.csv"), Table, header=8)
dat2 = CSV.read(datadir("exp_raw", "BSA", "BSA_Single_Shelf_CONV", "2023-09-12-13_REVO_AD_T.csv"), Table, header=1)

lyo_conv_full_pre = map(revo_columnrename, dat1)
time1 = range(0u"hr", step=1u"minute", length=length(lyo_conv_full_pre))
lyo_conv_full = Table(lyo_conv_full_pre; t=time1)

trim_PD = lyo_conv_full.phase .== 6
istart_lyo_pd = findfirst(trim_PD)
istart_lyo_pd = findlast(trim_PD)

lyo_conv_pd = Table(lyo_conv_full[trim_PD])
lyo_conv_pd.t .-= lyo_conv_pd.t[1]

# Plot temperature and pressure to make sure things look right
# plT = @df lyo_conv_pd plot(:t, :Tsh_sp)
# plp = @df lyo_conv_pd plot(:t, :pch_pir, ylim=(0, 300))
# @df lyo_conv_pd plot!(plp, :t,  :pch_cm)
# plot(plT, plp, layout=(2,1))

# Continue reading data
thm_conv_full = map(fiberoptic_columnrename, dat2)

# Conditions for identifying start and end of primary drying
istart_thm_pd = argmin(thm_conv_full.T2[3100:end]) + 3100
istart_thm_pd += 150 # Manual fudge
iend_thm_pd = length(thm_conv_full)

thm_conv_sub = Table(thm_conv_full[istart_thm_pd:iend_thm_pd])
thm_conv_pd = Table(thm_conv_sub; t = uconvert.(u"hr", thm_conv_sub.tstamp .- thm_conv_sub.tstamp[1]))
thm_conv_pd.t .-= thm_conv_pd.t[1]

# Identify second derivative condition
pch_pir_sm = savitzky_golay(lyo_conv_pd.pch_pir, 291, 3, deriv=0).y 
pch_pir_d2 = savitzky_golay(lyo_conv_pd.pch_pir, 291, 3, rate=1u"minute", deriv=2).y
t_end = lyo_conv_pd.t[argmax(pch_pir_d2[500:end])+499] - lyo_conv_pd.t[1] +0.0u"hr"

# Sanity check for second derivative
# @df lyo_conv_pd plot(:t, :pch_pir)
# @df lyo_conv_pd plot!(:t, pch_pir_sm)
# tendplot!(t_end)

# Identify temperature series for fitting
# @df thm_conv_pd exptfplot(:t, :T1, :T2, :T3, :T4)

fitdat_center = @df thm_conv_pd PrimaryDryFit(:t, (:T2[:t .< 25u"hr"], :T2[:t .< 29u"hr"]); t_end)
fitdat_edge = @df thm_conv_pd PrimaryDryFit(:t, :T4[:t .< 20u"hr"])

safesave(datadir("exp_pro", "BSA_conv_processed.jld2"), @strdict lyo_conv_pd lyo_conv_full thm_conv_pd thm_conv_full t_end fitdat_center fitdat_edge)
end

# ------ --------
# BSA RF

let
rfproc = CSV.read(datadir("exp_raw", "BSA", "BSA_Single_Shelf_RF+3hr", "2023-10-06-08_REVO_AND.csv"), Table, comment="#")
rftherm = CSV.read(datadir("exp_raw", "BSA", "BSA_Single_Shelf_RF+3hr", "Temp_Data.csv"), Table, comment="#")

lyo_rf_full_pre = map(revo_columnrename, rfproc)
time1 = range(0u"hr", step=1u"minute", length=length(lyo_rf_full_pre))
lyo_rf_full = Table(lyo_rf_full_pre; t=time1)

is_PD = (lyo_rf_full.phase .== 6 ) .& (lyo_rf_full.step .∈ [[1,2,3]])
lyo_rf_pd = lyo_rf_full[is_PD]
lyo_rf_pd.t .-= lyo_rf_pd.t[1]
# @df lyo_rf_pd plot(:t, :pch_pir)
# @df lyo_rf_pd plot(:t, :Tsh_sp)

# Identify second derivative condition
pch_pir_sm = savitzky_golay(lyo_rf_pd.pch_pir, 91, 3, deriv=0).y 
pch_pir_d2 = savitzky_golay(lyo_rf_pd.pch_pir, 91, 3, rate=1u"minute", deriv=2).y
t_end = lyo_rf_pd.t[argmax(pch_pir_d2[500:end])+499]  +0.0u"hr"

# Sanity check for second derivative
# @df lyo_conv_pd plot(:t, :pch_pir)
# @df lyo_conv_pd plot!(:t, pch_pir_sm)
# tendplot!(t_end)

# Start of cycle in earnest, for temperature data file
thm_rf_full_pre = map(fiberoptic_columnrename, rftherm)
thm_rf_full = Table(thm_rf_full_pre; t = uconvert.(u"hr", thm_rf_full_pre.tstamp))
# @df thm_rf_full exptfplot(:t, :T1, :T2, :T3, :T4)

istart_th = argmin(thm_rf_full.T4[3100:end]) + 3100
# tendplot!(thm_rf_full.t[istart_th])

thm_rf_pd = thm_rf_full[istart_th:end]
thm_rf_pd.t .-= thm_rf_pd.t[1]
# @df lyo_rf_pd plot(:t, :Tsh_sp, c=5)
# @df thm_rf_pd exptfplot!(:t, :T1, :T2, :T3, :T4)
# tendplot!(t_end)

# Identify data for fitting
fitdat_rf = @df thm_rf_pd PrimaryDryFit(:t, (:T2[:t .< 10u"hr"],
    :T3[:t .< 10u"hr"]), :T4, t_end)
fitdat_rf_alt = @df thm_rf_pd PrimaryDryFit(:t, (:T2[:t .< 10u"hr"],
    :T3[:t .< 10u"hr"]), :T3[findlast(:t.<t_end)] , t_end)

safesave(datadir("exp_pro", "BSA_rf_processed.jld2"), @strdict lyo_rf_pd lyo_rf_full thm_rf_full thm_rf_pd t_end fitdat_rf fitdat_rf_alt)
end

# --------------------------------
# M4 case, with variable power and MTSFlex sensors

let 
dat1 = CSV.read(datadir("exp_raw", "M4", "2023-02-08-09_LS_AD_2.csv"), Table, 
    dateformat="mm/dd/yyyy H:M",comment=";", stripwhitespace=true)
dat2 = CSV.read(datadir("exp_raw", "M4", "2023-02-08-09_LS_AD_1.csv"), Table, comment="#")
dat3 = CSV.read(datadir("exp_raw", "M4", "2023-02-08-09_LS_AD_3.csv"), Table, comment="#", stripwhitespace=true)

time2 = dat2.var"Elapsed [s]"

lyo_full_pre = map(lyostar_columnrename, dat1)
time1 = uconvert.(u"hr", lyo_full_pre.tstamp .- lyo_full_pre.tstamp[1])
lyo_full = Table(lyo_full_pre; t=time1)

# Conditions for identifying start and end of primary drying data
istart_lyo_pd = findfirst(lyo_full.pch_sp .== 100u"mTorr")
iend_lyo_pd = findlast(lyo_full.Tsh_sp .<= 10u"°C")

lyo_pd = Table(lyo_full[istart_lyo_pd:iend_lyo_pd])
lyo_pd.t .-= lyo_pd.t[1]

# Plot temperature and pressure to make sure things look right
plT = plot(lyo_pd.Tsh_sp)
plp = plot(lyo_pd.pch_pir, ylim=(0, 300))
plot!(plp, lyo_pd.pch_cm)
plot(plT, plp, layout=(2,1))

# Identify second derivative condition
pch_pir_sm = savitzky_golay(lyo_pd.pch_pir, 91, 3, deriv=0).y 
pch_pir_d2 = savitzky_golay(lyo_pd.pch_pir, 91, 3, rate=1u"minute", deriv=2).y
t_end = lyo_pd.t[argmax(pch_pir_d2[100:end])+99] - lyo_pd.t[1] +0.0u"hr"

# Sanity check for second derivative
@df lyo_pd plot(:t, :pch_pir)
@df lyo_pd plot!(:t, :pch_cm)
@df lyo_pd plot!(:t, pch_pir_sm)
tendplot!(t_end)

# Continue reading data
thm_full = map(fiberoptic_columnrename, dat2)

# Conditions for identifying start and end of primary drying
istart_thm_pd = 1
iend_thm_pd = length(thm_full)
istart_thm_pd = findlast(thm_full.T3 .< lyo_pd.Tsh_sp[1]) 
iend_thm_pd = findlast(thm_full.T3 .< 50u"°C")

thm_pd_sub = Table(thm_full[istart_thm_pd:iend_thm_pd])
thm_pd = Table(thm_pd_sub; t = uconvert.(u"hr", thm_pd_sub.tstamp .- thm_pd_sub.tstamp[1]))

# Plot temperatures for help in identifying type and useful region
@df thm_pd exptfplot(:t, :T3,)
@df thm_pd exptvwplot!(:t, :T1, :T4)
@df lyo_pd plot!(:t, :Tsh_sp, legend=:bottomright)

# Identify all the times when the RF is on and off
t_off = [0.0u"hr"]
t_on = [0.0u"hr"]
for _ in 1:11
    window = findall(t_on[end] .< thm_pd.t .< t_on[end]+45u"minute")
    # imax = findfirst((thm_pd.T3[window .+ 1] .- thm_pd.T3[window]) .< -0.1u"K") + window[1]
    imax = argmin(savitzky_golay(thm_pd.T3[window].|>u"K", 41, 3, deriv=2).y) + window[1]
    imin = argmin(thm_pd.T1[imax: imax+50]) + imax+1
    push!(t_off, thm_pd.t[imax])
    push!(t_on, thm_pd.t[imin])
end
@df thm_pd exptfplot(:t, :T3,)
@df thm_pd exptvwplot!(:t, :T1, :T4)
vline!(t_off, alpha=0.5)
vline!(t_on, alpha=0.5)

fitdat = @df thm_pd[begin:6:end] PrimaryDryFit(:t, :T3[:t.<t_off[5]], (:T1, :T4), t_end)
plot(fitdat)

# # Special: power has some modulation
tstops = sort(vcat(t_on, t_off))[2:end]
powers = [(t∈t_on ? (t < 5.8u"hr" ? 18u"W" : 12.5u"W")*0.54/17 : 0.0u"W" ) for t in tstops]
P_per_vial = ConstantInterpolation(powers, tstops, dir=:left, extrapolation_right=ExtrapolationType.Constant)
# P_per_vial_an = t -> (t < 5.8u"hr" ? 18u"W" : 12.5u"W")/17 * (t > t_on[findlast(t_off .<= t)])
# pow_full = map(powermeter_columnrename, dat3)
# istart_pow = findlast(pow_full.tstamp .<= 5u"s")
# pow_pd = pow_full[istart_pow:end]
# P_per_vial = LinearInterpolation(pow_pd.P .*0.54./17, pow_pd.tstamp, extrapolation=ExtrapolationType.Linear) # Correct for efficiency, # of vials

safesave(datadir("exp_pro", "M4_processed.jld2"), @strdict lyo_full lyo_pd thm_pd thm_full t_end fitdat P_per_vial)

end

# ------------------------------------
# SM1: conv

let 
dat1 = CSV.read(datadir("exp_raw", "SM1", "FiberOptic_27Apr2023.csv"), Table; comment=";")
dat2 = CSV.read(datadir("exp_raw", "SM1", "4_27_Con_M_S.csv"), Table, 
    comment=";", stripwhitespace=true, dateformat="mm/dd/yyyy H:M:S")

time1 = dat1.var"Elapsed [s]"

lyo_full_pre = map(lyostar_columnrename, dat2)
time2 = uconvert.(u"hr", lyo_full_pre.tstamp .- lyo_full_pre.tstamp[1])
lyo_conv_full = Table(lyo_full_pre; t=time2)

# Conditions for identifying start and end of primary drying data
istart_lyo_pd = findfirst(lyo_conv_full.pch_sp .== 70u"mTorr")
iend_lyo_pd = length(lyo_conv_full)

lyo_conv_pd = Table(lyo_conv_full[istart_lyo_pd:iend_lyo_pd])
lyo_conv_pd.t .-= lyo_conv_pd.t[1]

# Plot temperature and pressure to make sure things look right
# plT = plot(lyo_conv_pd.Tsh_sp)
# plp = plot(lyo_conv_pd.pch_pir, ylim=(0, 300))
# plot!(plp, lyo_conv_pd.pch_cm)
# plot(plT, plp, layout=(2,1))

# Continue reading data
thm_conv_full = map(fiberoptic_columnrename, dat1)

# Conditions for identifying start and end of primary drying
istart_thm_pd = argmin(thm_conv_full.T2[3100:end]) + 3099
iend_thm_pd = length(thm_conv_full)

thm_pd_sub = Table(thm_conv_full[istart_thm_pd:iend_thm_pd])
thm_conv_pd = Table(thm_pd_sub; t = uconvert.(u"hr", thm_pd_sub.tstamp .- thm_pd_sub.tstamp[1]))

# Identify end of drying
t_end = lyo_conv_pd.t[end] +0.0u"hr" # Just at end of data is basically the best point 

# Sanity check for drying end
# @df lyo_conv_pd plot(:t, :pch_pir)
# @df plot!(:t, pch_pir_sm)
# tendplot!(t_end)

# Plot temperatures for help in identifying type and useful region
# @df thm_pd exptfplot(:t, :T2, :T4)
# plot!(xlim=(0, 10))

iend_T2 = findlast(thm_conv_pd.t .< 8.5u"hr")
iend_T4 = findlast(thm_conv_pd.t .< 7u"hr")

# Condense down into the part that will be fit
fitdat_conv = @df thm_conv_pd PrimaryDryFit(:t, (:T2[begin:iend_T2], 
    :T4[begin:iend_T4]); t_end);

# Save data in its postprocessed form, to avoid repeating the above
safesave(datadir("exp_pro", "SM1_processed.jld2"), @strdict lyo_conv_pd thm_conv_pd thm_conv_full lyo_conv_full t_end fitdat_conv)
    
end

#### -------------------------
# SM2: RF
let
dat2 = CSV.read(datadir("exp_raw", "SM2", "5_7_RF1_M_S.csv"), Table, 
    comment=";", stripwhitespace=true, dateformat="mm/dd/yyyy H:M:S")
dat1 = CSV.read(datadir("exp_raw", "SM2", "FiberOptic_06May2023.csv"), Table)

time1 = dat1.var"Elapsed [s]"

lyo_full_pre = map(lyostar_columnrename, dat2)
time2 = uconvert.(u"hr", lyo_full_pre.tstamp .- lyo_full_pre.tstamp[1])
lyo_rf_full = Table(lyo_full_pre; t=time2)

istart_lyo_pd = findfirst(lyo_rf_full.pch_sp .== 70u"mTorr")
iend_lyo_pd = length(lyo_rf_full)-100 # Trim off last bit of data

lyo_rf_pd = Table(lyo_rf_full[istart_lyo_pd:iend_lyo_pd])
lyo_rf_pd.t .-= lyo_rf_pd.t[1]

# Plot temperature and pressure to make sure things look right
plT = @df lyo_rf_pd plot(:t, :Tsh_sp)
plp = plot(lyo_rf_pd.pch_pir, ylim=(0, 300))
plot!(plp, lyo_rf_pd.pch_cm)
plot(plT, plp, layout=(2,1))

# Identify end of drying
pch_pir_sm = savitzky_golay(lyo_rf_pd.pch_pir, 91, 3, deriv=0).y 
pch_pir_d2 = savitzky_golay(lyo_rf_pd.pch_pir, 91, 3, rate=1u"minute", deriv=2).y
t_end = lyo_rf_pd.t[argmax(pch_pir_d2[100:end])+99]  +0.0u"hr"

# Sanity check for second derivative
# @df lyo_rf_pd plot(:t, :pch_pir)
# @df lyo_rf_pd plot!(:t, pch_pir_sm)
# tendplot!(t_end)

# Start of cycle in earnest, for temperature data file
thm_rf_full_pre = map(fiberoptic_columnrename, dat1)
thm_rf_full = Table(thm_rf_full_pre; t = uconvert.(u"hr", thm_rf_full_pre.tstamp.-thm_rf_full_pre.tstamp[1]))
@df thm_rf_full exptfplot(:t, :T1, :T2, :T4)


# istart_thm_pd = argmin(thm_rf_full.T2[3100:end]) + 3099
# iend_thm_pd = length(thm_rf_full)

# thm_pd_sub = Table(thm_rf_full[istart_thm_pd:iend_thm_pd])
# thm_rf_pd = Table(thm_pd_sub; t = uconvert.(u"hr", thm_pd_sub.tstamp .- thm_pd_sub.tstamp[1]))

istart_th = @df thm_rf_full argmin(:T2[20u"hr" .< :t .< 25u"hr"])
istart_th += findlast(thm_rf_full.t .< 20u"hr")
iend_th = @df thm_rf_full findlast(:t .< t_end+:t[istart_th])

thm_rf_pd = thm_rf_full[istart_th:iend_th]
thm_rf_pd.t .-= thm_rf_pd.t[1]
# @df lyo_rf_pd plot(:t, :Tsh_sp, c=5)
# @df thm_rf_pd exptfplot!(:t, :T1, :T2, :T4)
# tendplot!(t_end)

#--------- Trim temperature data for fitting

fitdat_rf = @df thm_rf_pd PrimaryDryFit(:t, (:T2[:t.<3u"hr"], 
    :T4[:t.<3u"hr"]), :T1[begin:argmax(:T4)], t_end)
# plot(fitdat_rf)

safesave(datadir("exp_pro", "SM2_processed.jld2"), @strdict lyo_rf_pd lyo_rf_full thm_rf_full thm_rf_pd t_end fitdat_rf )
end
