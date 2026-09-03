# DFSEM Manual

Preparation: load **OpenFOAM v2012**

1. decomposePar
2. mpirun -np 8 pimpleFoam -parallel

The open, suburban, and urban case can be implemented by replacing the L, R, U-files in constant>boundaryData>inlet>0.
