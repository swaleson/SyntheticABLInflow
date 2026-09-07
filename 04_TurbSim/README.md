#TurbSim Manual
Preparation: Download and install the TurbSim executable (https://www.nlr.gov/wind/nwtc/turbsim) and adjust the TurbSim.inp file according to the adjustments for the open, suburban, and urban case as provided in the text files.. Beware of the data disclaimer (https://www.nlr.gov/disclaimer) and that the NREL was changed into NLR recently. The structure of the case is based on the DFSR example case provided by Melaku & Bitsuamlak (2021)

1.	Open a command window, navigate to the location and execute TurbSim to generate the flowfield.
2.	Split the flow field into OpenFOAM time steps by running xxx.m from Section 4 of this repository, ensure the data is stored in constant>boundaryData.
3.	decomposePar
4.	mpirun -np 8 pimpleFoam -parallel
