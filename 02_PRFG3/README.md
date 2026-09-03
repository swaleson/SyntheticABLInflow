# PRFG³ Manual

Preparation: After the .zip is downloaded and unpacked, the directories *CodeOF*, *inflowCorrection*, and *inflowGeneration* are copied into the base_case. Copy mesh into inflowCorrection directory. 

1. python inflowGeneration.py 
2. python inflowCorrection.py
3. wmake libso (with **OpenFOAM v2112** loaded)

Copy inflowCorrectionDict and inflowGenerationDict into system directory, copy synInflowLib.so into base_case

4. decomposePar
5. mpirun -np 8 pimpleFoam -parallel

The open, suburban, and urban case can be implemented by replacing the inflowGeneration.py file in the inflowGeneration directory. 
