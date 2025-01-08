# Flux RETAP: A REaction Target Prioritization Genome-Scale Modeling Technique for Selecting Genetic Targets

A function for selecting reaction targets that have changes that are correlated/anti-correlated with product production. 

- [Setup](#setup)
- [Instructions for use](#instructions-for-use)
- [Summary of files](#summary-of-files)
- [Reference](#reference)
- [License](#license)
- [Copyright Notice](#copyright-notice)


## Setup

Flux RETAP utilizes the cobra toolbox to perform simulations and has the following dependences:

* [cobra](https://opencobra.github.io/cobrapy/)
* scipy
* pandas
* numpy
* matplotlib

If the required packages are not installed, you can use the following line from the command prompt to install:
> pip install cobra scipy pandas numpy matplotlib

## Instructions for use

To run Flux RETAP, download the FluxRETAP.py function and move it into the directory containing your code. The module can be imported the following command:
>from FluxRETAP import FluxRETAP

Before simulation, it is necessary to import the [cobra package](https://opencobra.github.io/cobrapy/) and load the desired genome-scale model (see [cobra documentation](https://opencobra.github.io/cobrapy/) for more information)
It is also necessary to supply FluxRETAP with:
* product reaction name 
* carbon source reaction name
* biomass reaction name
* desired subsystems list (from the model)

A full tutorial demonstrating the parameters and use of FluxRETAP is provided [here](https://github.com/JBEI/FluxRETAP/blob/main/jupyter%20notebooks/FluxRETAP_Tutorial.ipynb) 


## Summary of files

* core - directory containing main code of FluxRETAP and FluxRETAP functions
* FluxRETAP_Tutorial.ipynb - notebook demonstrating how to use FluxRETAP
* FluxRETAP_Results.ipynb - results reported in the manuscript comparing flux retap to other GSM techniques and experimental identified genes
* Sensitivity.ipynb - sensitivity analysis of changing the parameters for heterologous indigoidine produciton in three different GSMs for _Pseudomonas putida_, _Sacchromyces cerevisiae_, and _Escherichia coli_

## Reference

This project is not yet published

## Copyright Notice

Flux REaction TArget Prioritization (Flux RETAP) Copyright (c) 2024, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.

## License

Flux REaction TArget Prioritization (Flux RETAP) Copyright (c) 2024, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches,
or upgrades to the features, functionality or performance of the source
code ("Enhancements") to anyone; however, if you choose to make your
Enhancements available either publicly, or directly to Lawrence Berkeley
National Laboratory, without imposing a separate written license agreement
for such Enhancements, then you hereby grant the following license: a
non-exclusive, royalty-free perpetual license to install, use, modify,
prepare derivative works, incorporate into other computer software,
distribute, and sublicense such enhancements or derivative works thereof,
in binary and source code form.


