---
# see https://joss.readthedocs.io/en/latest/example_paper.html
# Use plain text in the title; do not use backticks for code formatting.
title: 'LyoPronto.jl: A Julia package for pharmaceutical lyophilization simulation and parameter fitting'
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
  - name: Purdue University, United States
    ror: 02dqehb95
    index: 1
date: 13 August 2017
bibliography: paper.bib

---

# Summary

`LyoPronto.jl` is a Julia package for simulating and optimizing pharmaceutical lyophilization (freeze-drying) processes. It provides an efficient implementation of the most common model, and tools for working with new heat and mass transfer models and fitting model parameters to experimental data.
These parameters, measured at laboratory scale, are used to optimize process conditions and to carry out scale-up and tech transfer activities. 

# Statement of Need

Lyophilization, better known as freeze drying, is widely used in the pharmaceutical industry 
to stabilize drug products which have a short shelf life in the liquid state.
The number of lyophilized drugs approved by the US FDA has increased significantly in recent decades [@lyohub2025AnnualReport].
Best practice in lyophilization uses a model published in 1985 [@pikalUseLaboratoryData1985; @tchessalovBestPracticesGuidelines2022; @tchessalovPracticalAdviceScientific2023] for process development, scale-up, and tech transfer. This model is formally a differential-algebraic equation (DAE) system, although few practitioners explicitly recognize it as such. 

Additionally, due to the inherent physical limitations of conventional lyophilization,
a batch process which may last as much as a week for some products, 
a variety of novel processes are being actively researched which bear broad similarity.
These include (but are not limited to):

- microwave-assisted lyophililization (and similar processes sometimes called microwave vacuum drying) [@hardterMicrowaveAssistedFreezeDrying2023]
- thin-film freeze drying [@southardMechanisticModelingPrimary2026]
- single-vial continuous freeze drying [@schaalOptimizationContinuousSpinfreezedrying2025]
- spin freeze drying [@capozziBatchContinuousFreezeDrying2019]
- spray freeze drying [@borgessebastiaoAtmosphericSprayFreezeDrying2017]
- "lyobeads", millimeter-scale droplets frozen and dried, which are the subject of ongoing work at LyoHUB 

For process design and optimization, these novel processes will require new physical models, 
which in many cases will be broadly similar to that for conventional freeze drying.
`LyoPronto.jl` is designed both for practitioners of conventional lyophilization 
and to facilitate research in models for novel processes.

`LyoPronto.jl` was built in Julia to leverage its capabilities for scientific computing.
First, the solvers for differential-algebraic equations (DAEs) provided by `OrdinaryDiffEq.jl` [@rackauckas2017differentialequations; @steinebachConstructionRosenbrockWanner2023], particularly the `Rodas5P()` method, are a much more natural fit than plain ODE solvers for the typical model of lyophilization, which is a DAE. This allows for more efficient simulation of the process without requiring custom implementation of the DAE solver.
Second, automatic differentiation (AD), particularly the `ForwardDiff.jl` package
[@revelsForwardModeAutomaticDifferentiation2016],
allows for efficient and accurate computation of gradients of simulation outputs with respect to parameters, which is important for parameter fitting. 
Third, `Unitful.jl` checks for correct physical dimensionals and efficiently computes physical quantities with units, which is important for a field that frequently reports quantities such as $14.0\,\mathrm{cm^2\,Torr\,hr/g}$ (equal to $6.71 \times 10^5\, \mathrm{m/s}$)  for the mass transfer resistance $R_p$.
Fourth, `Plots.jl` [@tombreloffPlotsjl2026] provides a natural way to provide sensible default visualizations of experiment and simulation results, while allowing for significant customization.
Finally, Julia's performance characteristics allow for faster simulation and optimization than would be possible with interpreted languages like Python.


# State of the Field

Several other tools exist for solving the main literature model for conventional lyophilization:                                                     

- The term "LyoCalculator" is often used to refer to any of several Excel templates which implement the most common model for lyophilization. [One such template](https://pharmahub.org/groups/lyo/File:Lyocycle_design_and_transfer_template_Breckenridge_workshop--Pfizer_Original.xls) is still [distributed by LyoHUB](https://pharmahub.org/groups/lyo/tools) at Purdue University. By using Excel, these tools are made accessible to a wide audience (particularly to non-engineers), but it is painstaking to verify correctness if a new model is developed (and impossible to extend for the digital twins of the future).
- LyoPRONTO (note the all-caps), a Python package, was previously developed by other researchers at LyoHUB [@shivkumarLyoPRONTOOpenSourceLyophilization2019]. This package has received some recent updates, but employs very simplistic numerical approaches to simulation and parameter fitting (which motivated the development of `LyoPronto.jl`).
- A web interface to LyoPRONTO was created in 2019, and is currently available at [https://lyopronto.geddes.rcac.purdue.edu](https://lyopronto.geddes.rcac.purdue.edu). This interface exposes all the main functionality available in LyoPRONTO at its release, but does not allow deep inspection of results (such as intermediate calculations) or extension of the code base. All that could be found of its web frontend's source is minified JavaScript, which makes maintenance difficult.
- A revised web interface to LyoPRONTO which adds some functionality was recently made available at [https://lyopronto2.geddes.rcac.purdue.edu](https://lyopronto2.geddes.rcac.purdue.edu) (note the 2). This interface is written in Python using Streamlit, which should make it more maintainable for scientists. This interface exposes some functionality from `LyoPronto.jl`.
- Another software tool was announced with a journal article by @bognerSoftwareToolLyophilization2021. This tool is not open source, and the journal article does not give a name for the software or explain how to obtain it, short of contacting the authors directly.
- Several large pharmaceutical companies maintain their own internal software tools for lyophilization process design and analysis, but these are not publicly available.

# Software Design

`LyoPronto.jl` is designed to be flexible and extensible, particularly along three axes:

1. Adding new models for the lyophilization process or for similar processes, such as microwave-assisted lyophilization. 
2. In fitting, comparing different types or amounts of experimental data to simulation, such as such as product temperatures, vial wall temperatures, and drying times. Some of these are not always available in a single dataset, but including as many as are available will improve the quality of the fit.
3. Allowing different sets of model parameters to be fit to without redefining the entire model and experimental datasets. 

To achieve this flexibility,
the loss function for parameter fitting is defined
to take a plain vector of numbers as input, along with a 
transform, a struct of base model parameters, and the set of experimental data.

The transform is used to map the input vector 
to a named tuple of model parameters with physical units
using `TransformVariables.jl`.
The named parameters are then used to create a copy of the full set of simulation parameters with the new values, leaving the others unchanged. 
This struct of parameters is then passed to the simulation function, which returns a 
completed simulation. 
The residuals between the simulation and experiment
are squared, reduced to a dimensionless form by user-defined weights 
(e.g., divide temperature loss by $1 \mathrm{K}^2$ and divide the drying time loss by $1 \mathrm{hr}^2$) 
and then summed to yield the loss function value as a plain scalar.

Thus, to change the set of parameters to be fit (or the range of their allowed values), it is enough to change the transform; to use different experimental data, simply pass a different dataset; to change the model, provide a new struct of parameters (to match a new model structure). The loss function itself does not need to be redefined, and the simulation code does not need to be modified.

`LyoPronto.jl`'s online documentation provides a tutorial for adding new models and for adding new types of experimental data. 
There is also a demonstration of how to extend the transform framework  
to fit multiple experiments at once with some fit parameters in common and others separate for each experiment, 
such as if three experiments are expected to share heat transfer behavior
but have different mass transfer behavior.


# Research Impact Statement

Our prior work on microwave-assisted lyophilization
[@alexeenkoRandomizedfieldMicrowaveassistedPharmaceutical2025] 
used an implementation in Python of a model for microwave-assisted lyophilization, 
but least-squares fitting with that framework was difficult so all fit parameters were manually tuned by visual comparison between experiment and simulation.

With `LyoPronto.jl`, we can fit an improved microwave model to experimental data using a rigorous least-squares criterion.
In [@wheelerMechanisticModelsExperimental2026], 
we usually had access to temperature measurement of the product in at least one vial,
but only sometimes had corresponding measurements of the vial wall.
We were able to determine that product temperature at the end of drying tends to match the final vial wall temperature, 
adding another type of experimental data to the the fitting process.
With `LyoPronto.jl`, we could apply all of the available data to each experimental fit, 
using a rigorous least-squares criterion for fitting,
without redefining the loss function or modifying the simulation code.

Likewise, in [@wheelerMechanisticModelsExperimental2026], we usually fit only the microwave-related parameters of the model
to the microwave experimental data,
with mass transfer parameters fit to conventional lyophilization data.
However, in some cases, other evidence suggested that the mass transfer parameters
would be different in the microwave experiments, and `LyoPronto.jl` could fit both microwave and mass transfer parameters to a single experiment by simply using a different transform in parameter fitting.

In work on a population-balance model for pore size evolution during lyophilization [@wheelerDirectModelingPore2026], we were able to use `LyoPronto.jl` to fit a new model with new kinds of parameters to experimental data, including direct comparisons to the conventional model.

Ongoing work on a model for the lyophilization of "lyobeads" 
(millimeter-scale beads of frozen solution)
has already begun, extending `LyoPronto.jl` with a new model 
and fitting that model to experimental data.

Although not the primary focus of this software, the equipment capability model constructed by 
@kazarinCompactModelLyophilizer2022 
(originally provided as a MATLAB script in the supplementary material) 
is also implemented in `LyoPronto.jl` as an `ECCURT` module,
which is part of the backend for the improved web tool at [lyopronto2.geddes.rcac.purdue.edu](http://lyopronto2.geddes.rcac.purdue.edu).
This web interface is provided to the public by LyoHUB, and is used by researchers doing process design and optimization for conventional lyophilization. 
`LyoPronto.jl` is actively being incorporated to provide other functionality for that web interface, especially for parameter fitting.

# AI usage disclosure

For most of its history, `LyoPronto.jl` has been developed without the use of generative AI tools.
Recent development has used autocomplete from GitHub Copilot, 
and some PRs were authored
with assistance from Qwen3.8 running on university servers; 
all code was reviewed and revised by the authors. 
This manuscript was written without the use of generative AI tools.

# Acknowledgements

Part of the scientific work which motivated development of `LyoPronto.jl` was funded by the 
National Institute for Innovation in Manufacturing Biopharmaceuticals (NIIMBL) under project PC4.1-307.

My PhD advisors, Prof. Alina Alexeenko and Prof. Vivek Narsimhan, provided scientific guidance and support, although they were not directly involved in the development of `LyoPronto.jl`.

# References