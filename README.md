# InCens
Information transfer in Central Dogma systems

This package contains python3 scripts and jupyter notebooks to compute the ideal and the protein-level channel capacity of central dogma reaction networks.
Code for Gillespie (kinetic Monte Carlo) simulation of central dogma reactions is in the scripts subdirectory, CentralDogmaSimulator.py. The code for computing 
the maximum information transfer through central dogma reactions, the ideal channel capacity, is also in the scripts subdirectory, IdealIntegration_analytical.py.

The central dogma rate constants for the four species are in scripts/simulation_data.py. To perform a Gillespie simulation the command is

$ python3 CentralDogmaSimulator.py -d DATANAME -n NUMBER OF SAMPLES

The output of the simulation, the protein expression values for each integration time, is stored in the subdirectory, simulation_results/DATANAME_samples/.
