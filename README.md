# Sucrose_breath_test_model

This repository contains data and example code for a mechanistic model of the 13C sucrose breath test (13C-SBT).
It is associated with the paper "Connecting 13C-sucrose breath test curve dynamics to underlying metabolic processes: development of a model-based 13C-sucrose breath test diagnostic for gut function disorders characterized by a loss of sucrase-isomaltase enzymatic activity" by Andrew F. Brouwer1, Gwenyth O. Lee, Hannah Van Wyk, Robert J. Schillinger, Christine A. Edwards, and Douglas J. Morrison

This repository includes two data sets
*DATA_SFG: Data from a collection of three 13C-SBT experiments. In the first experiment (Experiment S),  the sucrose tracer consisted of a highly and uniformly enriched fructose moiety (U-13C fructose) with an unlabeled glucose moiety; in the second experiment (Experiment G), the sucrose tracer consisted of a highly and uniformly enriched glucose moiety (U-13C glucose) with an unlabeled fructose moiety; in the last experiment (Experiment F), the sucrose tracer was highly and uniformly enriched (U-13C sucrose). The columns in the data indicate the participant, the experiment, the time, and the percent dose recovery rate (PDRr) in 1/hr.
*DATA_MLE: Data from a collection of three 13C-SBT experiments. In each of the three experiments, the sucrose tracer was U-13C sucrose. In addition to the sucrose tracer, participants were given a dose of mulberry leaf extract (MLE), Reducose®, standardized to contain 5% 1-Deoxynojirimycin, which is an α-glucosidase inhibitor. The mulberry leaf extract doses were 0 mg, 100 mg, and 750 mg for the three experiments. The columns in the data indicate the participant, the experiment, the time, and the percent dose recovery rate (PDRr) in 1/hr.

As an illustration, we provide code to fit the data from Experiment S using the mechanistic model and 
*Example code. Includes model, simulation, optimization, and plotting code.
*Results S. Participant number, number of data points (n), best fit negative log-likelihood (NLL) and Schwarz Information Criterion (SIC), model parameters on the transformed scale used for optimization (par[1], par[2], par[3]), and model parameters (rho, pi, kappa).
*Plots_S. Includes data (points) and best-fit model (line).

