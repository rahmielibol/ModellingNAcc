# Modelling Nucleus Accumbens
A Computational Model from Single Cell to Circuit Level

A computational model is established for nucleus accumbens in Python, using BRIAN2 library. The results are verified the validity of the model by showing the consistency of simulation results with the empirical data. So, to run the code, one should download BRIAN2 (https://briansimulator.org/)

First, single cells are considered and their model behavior, then synaptic currents are considered. The results are given with single neuron behavior (code lines 1048-1120), synaptic currents (code lines 1124-1210), raster plots (code lines 1215-1327) and local field potentials (code lines 1520-1550). There are two scenarios, scenario 0 is for testing whether all is working OK, and scenario 2 stimulus is given and the overall behavior of the nucleus accumbens is observed.

As there are random variables, the results of two runs would not be exactly same, but very similar and behavior would be consistent. All the results are written in files.

 Further detail is written as comment lines in the code.
