---
# see https://joss.readthedocs.io/en/latest/example_paper.html
# Use plain text in the title; do not use backticks for code formatting.
title: 'LyoPronto.jl: A Julia package for pharmaceutical lyophilization'
tags:
  - Julia
  - pharmaceutics
  - lyophilization
  - parameter fitting
  - 
authors:
  - name: Isaac S. Wheeler
    orcid: 0000-0002-9717-073X
    corresponding: true
    affiliation: 1 
affiliations:
 - name: Davidson School of Chemical Engineering, Purdue University, United States
   index: 1
date: 13 August 2017
bibliography: paper.bib

---

# Summary

`LyoPronto.jl` is a Julia package for simulating and optimizing pharmaceutical lyophilization (freeze-drying) processes. It provides tools for modeling the heat and mass transfer during lyophilization and fitting model parameters to experimental data.
These parameters, measured at laboratory scale, are used to optimize process conditions and to carry out scale-up and tech transfer activities. 

# Statement of Need

Best practice in lyophilization uses a model published in 1985 [@pikal1985; @tchessalov2022; @tchessalov2023] for process development, scale-up, and tech transfer. This model is formally a differential-algebraic equation (DAE) system, although few practitioners explicitly recognize it as such. 

Several other tools exist for solving this model lyophilization:                                                     
- The term "LyoCalculator" is often used to refer to any of several Excel templates which implement the most common model for lyophilization. [One such template](paper-resources/Lyocycle_design_and_transfer_template_Breckenridge_workshop--Pfizer_Original.xls) is still [distributed by LyoHUB](https://pharmahub.org/groups/lyo/tools) at Purdue University. By using Excel, these tools are made accessible to a wide audience (particularly to non-engineers), but it is painstaking to verify correctness if a new model is developed (and impossible to extend for the digital twins of the future).
- A software tool was announced with a paper in 2021 by @bogner2021. This tool is not open source, and the journal article does not make clear a title of the software or how to obtain it, short of contacting the authors directly.
- LyoPRONTO (note the all-caps), a Python package, was previously developed by other researchers at LyoHUB [@shivkumar2019]. This package is still partially maintained, but employs very simplistic numerical approaches to simulation and optimization (which motivated the development of `LyoPronto.jl`).
- A web interface to LyoPRONTO was hosted at launch, and currently is made available at [https://lyopronto.geddes.rcac.purdue.edu](https://lyopronto.geddes.rcac.purdue.edu). This interface exposes all the main functionality available in LyoPRONTO at its release, but does not allow deep inspection of results (such as intermediate calculations) or extension of the code base. All that could be found of its web interface source is minified JavaScript, which makes maintenance difficult.
- A revised web interface to LyoPRONTO which adds some functionality was recently made available at [https://lyopronto2.geddes.rcac.purdue.edu](https://lyopronto2.geddes.rcac.purdue.edu) (note the 2). This interface is written in Python using Streamlit, which should make it more maintainable for scientists. This interface also exposes some functionality from `LyoPronto.jl`.
- Several pharmaceutical companies maintain internal software tools for lyophilization, but these are not publicly available.
                                                                                                                                        
`LyoPronto.jl` was built in Julia to leverage its capabilities for scientific computing.
First, the solvers for differential-algebraic equations (DAEs) provided by `OrdinaryDiffEq.jl` [@rackauckas2017differentialequations; @steinebachConstructionRosenbrockWanner2023], particularly the `Rodas5P()` method, are a much more natural fit than pure ODE solvers for the typical model of lyophilization, which is a DAE. This allows for more efficient simulation of the process without requiring custom implementation of the DAE solver.
Second, automatic differentiation (AD), particularly the `ForwardDiff.jl` package [@Revels2016ForwardDiff], allows for efficient and accurate computation of gradients, which is important for parameter fitting. 
Third, `Unitful.jl` checks for correct physical dimensionals and efficiently computes physical quantities with units, which is important for a field that frequently reports quantities such as $14.0 \mathrm{cm^2\,Torr\,hr/g}$ (equal to $6.71 \times 10^5 \mathrm{m/s}$)  for the mass transfer resistance $R_p$.
Fourth, `Plots.jl` [@tombreloffPlotsjl2026] provides a natural way to provide sensible default visualizations of experiment and simulation results, while allowing for significant customization.
Finally, Julia's performance characteristics allow for faster simulation and optimization than would be possible with interpreted languages like Python.

## New Models and New Parameter Fitting

`LyoPronto.jl` is designed to be flexible and extensible, 
including for new models of processes which are similar but not identical to conventional lyophilization. 
For example, in microwave-assisted lyophilization, 
new heat transfer mechanisms are present, 
and the model must be modified to account for these.
`LyoPronto.jl` natively provides one such model for microwave-assisted lyophilization [@wheelerMechanisticModelsExperimental2026],
which is not available in any of the other tools mentioned above.
In developing this type of model, 
it is often desirable to fit varying sets of model parameters to available experimental data. 
`LyoPronto.jl` allows the user to specify which parameters to fit, 
and which to hold constant, 
in a straightforward manner by using `TransformVariables.jl`. 
This is a feature that is not available in the original LyoPRONTO package or any of the other tools mentioned above.

This extensibility is exploited in other work, such as a population-balance model which
accounts for dynamic changes in the pore size distribution of the frozen cake during lyophilization [@wheelerDirectModelingPore2026], 
and ongoing work developing a model for "lyobeads" (millimeter-scale beads of frozen solution which are lyophilized outside a vial).

# AI usage disclosure

For most of its history, `LyoPronto.jl` has been developed without the use of generative AI tools.
Recent development has used autocomplete from GitHub Copilot, and some PRs were authored
with assistance from Qwen3.8 running on university servers; all code was reviewed and edited by the authors. 
This manuscript was written without the use of generative AI tools.

# Acknowledgements

Part of the scientific work which motivated development of `LyoPronto.jl` was funded by the 
National Institute for Innovation in Manufacturing Biopharmaceuticals (NIIMBL) under project PC4.1-307.

My PhD advisors, Prof. Alina Alexeenko and Prof. Vivek Narsimhan, provided scientific guidance and support, although they were not directly involved in the development of `LyoPronto.jl`.

# References