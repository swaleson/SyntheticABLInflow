# Mann Manual
Preparation: Download the method (https://github.com/fengguoFUAS/4D-Mann-Turbulence-Generator) and adjust the files MannTurbConfig.m and MannTurbFieldConfig.m as required. The open, suburban, and urban case can be implemented by the provided files. The structure of the case is based on the DFSR example case provided by Melaku & Bitsuamlak (2021)

1.	Generate the flowfield by running RunExampleAndVizualize.m which starts the executable
2.	Split the flow field into OpenFOAM time steps by running xxx.m from Section 4 of this repository, ensure the data is stored in constant>boundaryData.
3.	decomposePar
4.	mpirun -np 8 pimpleFoam -parallel

MIT License

Copyright (c) 2022 fengguoFUAS

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*
