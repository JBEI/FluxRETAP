# Flux RETAP: A REaction Target Prioritization Genome-Scale Modeling Technique for Selecting Genetic Targets

A function for selecting reaction targets that have changes that are correlated/anti-correlated with product production. 

- [Setup](#setup)
- [Instructions for use](#instructions-for-use)
- [Summary of files](#summary-of-files)
- [Reference](#reference)
- [License](#license)


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

## License

This code is currently under construction, and a license has not been chosen
