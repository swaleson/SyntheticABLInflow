# DFSR Manual

After the installation of DFSR, die **OpenFOAM v8** case is decomposed. DFSRTurb is run to generate the inflow data for each time step. Finally, the simulation can be started. The structure of the case is based on the DFSR example case provided by Melaku & Bitsuamlak (2021).

1. decomposePar
2. mpirun -np 8 DFSRTurb -parallel
3. mpirun -np 8 pimpleFoam -parallel

The open, suburban, and urban case can be implemented by adjusting the *profile*-file in constant>boundaryData>windProfile, by taking the respective file from the directory.
