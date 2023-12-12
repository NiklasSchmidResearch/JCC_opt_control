This repository contains code used to generate examples for the paper "Computing Optimal Joint Chance Constrained Control Policies" by Niklas Schmid, Marta Fochesato, Sarah H.Q. Li, Tobias Sutter and John Lygeros.

There are two simulations on a quadcopter example (A_main_quadcopter_example.m) and one simulation for a fishery management example (A_main_fishSimulation.m). The two quadcopter simulations can be adjusted via quadcopter_exmaple_settings.m. The settings for the fishery management simulation can be found in the beginning of the file A_main_fishSimulation.m. 

The files "filterMaskAndComputeKernel" are used to generate the transition kernel of the respective stochastic systems and a mask indicating safe and unsafe states.

The helper functions are used to run DP recursions to compute the cheapest, safest and "safest-under-Boole's-Inequality", i.e., the safest policy when safety is computed under Boole's Inequality, policies.

To run the code, one needs to also download the arrow3 files by Tom Davis available online at https://de.mathworks.com/matlabcentral/fileexchange/14056-arrow3
since they are used for one of the plots.
