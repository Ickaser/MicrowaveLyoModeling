# Preprocess experiments to get experimental data to choose usable data
let include(scriptsdir("postprocess_sugars.jl")) end

# Fit LC model to experiments
let include(scriptsdir("tuning_M1_RF.jl")) end
let include(scriptsdir("tuning_M3_RF.jl")) end
let include(scriptsdir("tuning_M4.jl")) end
let include(scriptsdir("tuning_SM_KvRpRF.jl")) end
# Fit LC model to literature experiments
let include(scriptsdir("bhambhani2021_modelcomparison.jl")) end
let include(scriptsdir("gitter2019_modelcomparison.jl")) end
# Generate figures which compare across cases
let include(scriptsdir("comparing_RF_cases_model.jl")) end
