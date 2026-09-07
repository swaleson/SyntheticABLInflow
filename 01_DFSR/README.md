# DFSR Manual
After the installation of DFSR, die OpenFOAM v8 case is decomposed. DFSRTurb is run to generate the inflow data for each time step. Finally, the simulation can be started. The structure of the case is based on the DFSR example case provided by Melaku & Bitsuamlak (2021).

1.	decomposePar
2.	mpirun -np 8 DFSRTurb -parallel
3.	mpirun -np 8 pimpleFoam -parallel

The open, suburban, and urban case can be implemented by adjusting the profile-file in constant>boundaryData>windProfile, by taking the respective file from the directory. The original files are from Melaku, A.F. and Bitsuamlak, G.T., 2021. A divergence-free inflow turbulence generator using spectral representation method for large-eddy simulation of ABL flows. Journal of Wind Engineering and Industrial Aerodynamics, 212, p.104580 and include the adapted values for this comparison.

The following license applies to DFSR: 

Copyright (C) 2021

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or(at your option) any later version.

Over and above the legal restrictions imposed by this license, if you use this software for an academic publication then you are expected to provide proper attribution. This can be to this code directly,

A.F. Melaku and G.T. Bitsuamlak, DFSR, v1.0 (2021). https://zenodo.org/badge/latestdoi/300517389.

or the original article

A.F. Melaku and G.T. Bitsuamlak. https://doi.org/10.1016/j.jweia.2021.104580.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
