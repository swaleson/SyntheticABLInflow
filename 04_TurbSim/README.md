# TurbSim Manual

Preparation: Download and install the TurbSim executable and adjust the TurbSim.inp file. The open, suburban, and urban case can be implemented by the provided files.
Beware that the NREL was changed into NLR recently. 

1. Open a command window, navigate to the location and execute TurbSim to generate the flowfield.
2. Split the flow field into OpenFOAM time steps by running xxx.m from Section 4 of this repository, ensure the data is stored in constant>boundaryData.
3. decomposePar
4. mpirun -np 8 pimpleFoam -parallel
