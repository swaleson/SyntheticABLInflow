# SyntheticABLInflow
Benchmark OpenFOAM cases for synthetic inflow generation methods (DFSR, PRFG³, TurbSim, Mann, DFSEM) for urban/suburban/open terrain.

# Associated publication
This repository contains the OpenFOAM cases, Matlab scripts, and data underlying the comparison in:

*Waleson S., Geleta T.N., Hartz C., Bitsuamlak G. 2026. Turbulent inflow generation for atmospheric boundary layer flows (II): Comparative study of synthetic methods. Journal of Wind Engineering an Industrial Aerodynamics. Volume 278, Article 106585, https://doi.org/10.1016/j.jweia.2026.106585*

It is a complementary part of the research output and no stand-alone description of the methods and methodology. 
If you use material from this repository, please cite the publication above. 

# Overview 
1. References to the synthetic inflow generation method sources (external)
2. OpenFOAM base cases (suburban)
3. input and reference data used in the comparison for the open, suburban, and urban case
4. Matlab script to generate inflow data from tensors

# 01. Synthetic inflow generation method sources
The third-party synthetic inflow generation methods are only provided as links. The external material remains subject to the terms and conditions of the respective providers.

**DFSR** 

Provided through https://github.com/abiyfantaye/DFSR with https://doi.org/10.5281/zenodo.5015100 and presented in Melaku, A.F. and Bitsuamlak, G.T., 2021. A divergence-free inflow turbulence generator using spectral representation method for large-eddy simulation of ABL flows. Journal of Wind Engineering and Industrial Aerodynamics, 212, p.104580. 

**PRFG³**

Provided through https://site.unibo.it/cwe-lamc/en/downloads/syninflow and presented in Patruno L., Ricci M. 2018. A systematic approach to the generation of synthetic turbulence using spectral method. Computer Methods in Applied Mechanics and Engineering 340, 881-904

**Mann**

Provided with adjusted eddy life-time and multiple fields of which only one is used through https://github.com/MSCA-LIKE/4D-Mann-Turbulence-Generator. Guo, F., Mann, J., 2022. Initial release of the 4d mann turbulence generator: Mannturb4d (version v1.0). http://dx.doi.org/10.5281/zenodo.6223785. Original method in Mann, J., 1994. The spatial structure of neutral atmospheric surface-layer turbulence. J. Fluid Mech. 273, 141–168.

**TurbSim** 

Provided through https://www.nlr.gov/wind/nwtc/turbsim. Beware that the NREL was renamed to NLR. 

**DFSEM**

*turbulentDFSEMInlet* for version 2012. e.g. https://doc.openfoam.com/2306/tools/processing/boundary-conditions/rtm/derived/inlet/turbulentDFSEMInlet/

# 02. OpenFOAM base case (suburban)
The OpenFOAM cases are provided for the suburban case and the five synthetic inflow generation methods in the respective directory. 

# 03. Input open / suburban / urban
The input data from the comparison targeting the ESDU open, suburban, and urban case are included in the directory. 

# 04. Matlab Tensor -> Inlet BC
The matlab script to deconstruct the tensors from Mann and TurbSim as inflow boundary condition per time step are included in a separate directory. 


