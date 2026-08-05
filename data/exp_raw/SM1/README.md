# Case SM1

This folder contains raw experimental data for case SM1 from the lyophilization experiments discussed in the article. 
This experiment had no microwave power, and is intended to provide a direct comparison to case SM2 which was microwave-assisted.

`4_27_Con_M_S.csv` provides the lyophilizer's process data. In keeping the same experimental procedure as the microwave-assisted cycles, the product temperature is recorded separately in `FiberOptic_27Apr2023.csv`.

To see how this data was postprocessed before model fitting, see `scripts/postprocess_sugars.jl`.

## Figure 
The left panel of following figure, reported in the article as Figure 14, shows the primary drying portion of this experiment and a model fit.

![figure showing experimental fit](SM1-2_combine_T_q.svg)