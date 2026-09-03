# Mann Manual

Preparation: Download the method and adjust the files MannTurbConfig.m and MannTurbFieldConfig.m as required. The open, suburban, and urban case can be implemented by the provided files.

1. Generate the flowfield by running RunExampleAndVizualize.m which starts the executable
2. Split the flow field into OpenFOAM time steps by running xxx.m from Section 4 of this repository, ensure the data is stored in constant>boundaryData.
3. decomposePar
4. mpirun -np 8 pimpleFoam -parallel
