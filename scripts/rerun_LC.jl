# All the below run inside `let` blocks to avoid variable name collisions
# between very-similar scripts.

# Preprocess experiments to get experimental data to choose usable data
let; include(scriptsdir("postprocess_sugars.jl")); end

# Compare LC variants
let; include(scriptsdir("M1_model_variants.jl")) end
# Fit LC model to experiments
let; include(scriptsdir("modelcomparison_M1.jl")) end
let; include(scriptsdir("modelcomparison_M3.jl")) end
let; include(scriptsdir("modelcomparison_M4.jl")) end
let; include(scriptsdir("modelcomparison_SM.jl")) end
# Fit LC model to literature experiments
let; include(scriptsdir("modelcomparison_bhambhani2021.jl")) end
let; include(scriptsdir("modelcomparison_gitter2019.jl")) end
# Generate figures which compare across cases
let; include(scriptsdir("comparing_RF_cases_model.jl")) end
