# DFSR Manual

After the installation of DFSR, die **OpenFOAM v8** case is decomposed. DFSRTurb is run to generate the inflow data for each time step. Finally, the simulation can be started. 

1. decomposePar
2. mpirun -np 8 DFSRTurb -parallel
3. mpirun -np 8 pimpleFoam -parallel

The open, suburban, and urban case can be implemented by adjusting the *profile*-file in constant>boundaryData>windProfile. 
