# DFSEM Manual
Preparation: load OpenFOAM v2012. The structure of the case is based on the DFSR example case provided by Melaku & Bitsuamlak (2021)

1.	decomposePar
2.	mpirun -np 8 pimpleFoam -parallel
   
The open, suburban, and urban case can be implemented by replacing the L, R, U-files in constant>boundaryData>inlet>0.
