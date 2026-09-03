# SyntheticABLInflow
Benchmark OpenFOAM cases for synthetic inflow generation methods (DFSR, PRFG³, TurbSim, Mann, DFSEM) for urban/suburban/open terrain.

# Associated publication
This repository contains the OpenFOAM cases, Matlab scripts, and data underlying the comparison in:

*Waleson S., Geleta T.N., Hartz C., Bitsuamlak G. 2026. Turbulent inflow generation for atmospheric boundary layer flows (II): Comparative study of synthetic methods. Journal of Wind Engineering an Industrial Aerodynamics. Volume 278, Article 106585, https://doi.org/10.1016/j.jweia.2026.106585*

It is a complementary part of the research output and no stand-alone description of the methods and methodology. 
If you use material from this repository, please cite the publication above. 

# Overview 
1. References to the synthetic inflow generation method sources (external)
2. OpenFOAM base case
3. input and reference data used in the comparison for the open, suburban, and urban case
4. Matlab script to generate inflow data from tensors

# 01. Synthetic inflow generation method sources
The third-party synthetic inflow generation methods are only provided as links. The external material remains subject to the terms and conditions of the respective providers.

| Inflow method | Reference | Availability |

| DFSR |  Melaku, A.F. and Bitsuamlak, G.T., 2021. A divergence-free inflow turbulence generator using spectral representation method for large-eddy simulation of ABL flows. Journal of Wind Engineering and Industrial Aerodynamics, 212, p.104580 | https://github.com/abiyfantaye/DFSR |

| PRFG³ | https://site.unibo.it/cwe-lamc/en/downloads/syninflow |

| TurbSim | https://www.nlr.gov/wind/nwtc/turbsim |

| Mann | https://github.com/MSCA-LIKE/4D-Mann-Turbulence-Generator | 
