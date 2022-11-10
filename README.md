# phybreak_multiple_introductions
Phybreak version and additional code used for the analyses in article of multiple introductions

The phybreak_1.0.0.tar.gz includes the package version used.

rm_mcmc.R is used for the simulation and mcmc run of outbreaks. 
Arguments for this function are the id number, saving directory, number of introductions in the simulation, outbreak size and the coalescent rate parameter of the history

simulation_analysis.R contains the code for the analysis of the mcmc runs. 
This code is used to produce figure 2 of the manuscript.

mink_analysis.R contains the code for preprocessing and mcmc run of the mink data.
This file also includes code for plots used in figure 3 and figures s5 and s6.

infectivity_function.R contains the code for the user-defined generation time distribution used with the mink farm analysis.
