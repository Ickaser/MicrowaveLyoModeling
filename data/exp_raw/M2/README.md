# Case M2

This folder contains raw experimental data for case M2 from the microwave-assisted lyophilization experiments discussed in the journal article. 

`RF_Petr_5_M_TEMP.csv` contains the process data recorded by the lyophilizer; since this was a microwave cycle, traditional thermocouples could not be used, and temperatures were logged on a separate computer.
`TemperatureLog_021023.csv` contains the fiber optic probe measurements of temperature, recorded from roughly the same time the cycle data logging began, although not exactly. 

To see how this data was postprocessed before model fitting, see `scripts/postprocess_sugars.jl`.

## Figure 
The following figure, reported in the article as Figure 18 (in the appendix), shows the primary drying portion of this experiment and a model fit.

![figure showing experimental fit](M2_T_q_combine.svg)