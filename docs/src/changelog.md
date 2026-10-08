# Changelog

## v0.5.0 - 2026-10-07

### ️ Breaking Changes

#### Reexports
- Refactored Re-exports: The package has moved away from heavy use of `@reexport`. Only the following symbols from other packages are still exported by `LyoPronto`; others (such as from `ConstructionBase`, `DiffEqCallbacks`, and `OrdinaryDiffEqNonlinearSolve`  must now be explicitly imported.
    - `ODEProblem` and `solve` from `OrdinaryDiffEqRosenbrock` remain explicitly exported.
        - Other `OrdinaryDiffEqRosenbrock` symbols are no longer re-exported, notably algorithms like `Rodas3()`.
        - The recommended algorithm for LyoPronto's ODE and DAE systems is now `Rodas5P(AutoForwardDiff(chunksize=2))`, which is made public as `LyoPronto.odealg_chunk2`.
    - From Unitful, the `u""` macro, `ustrip`, `uconvert`, and `NoUnits` remain explicitly exported.


#### Experimental Data & Fitting
- Experimental Data Handling: The `PrimaryDryFit` structure is now deprecated in favor of a new, modular hierarchy based on `ExpFitData` and `AbstractExpDatum`. 
    - Users should migrate from `PrimaryDryFit` to `ExpFitData` (e.g., using `TfData`, `TvwSeriesData`, etc.).
    - For help with migration, see the documentation for `ExpFitData` in the package API, refer to the `AbstractExpDatum` docstring, and consult the examples in the documentation.
- Parameter Objects: Model parameters are now documented as handled via structured, typed objects. Tuple-of-tuples data structures are only used for constructing the parameter objects.
- Loss and Residual Functions: the functions `obj_expT`, `err_expT`, and `err_expT!` are removed and replaced by `obj_exp`, `err_exp`, and `err_exp!`, which are more general and can handle multiple experimental data types. 
- Residual Weighting: The API for residual weighting has been updated to be more consistent across different experimental data types. The functions `loss_weighting` and `residual_weighting` are now used to specify weighting schemes for residuals in fitting, and these are passed to `obj_exp` and `err_exp!` functions. Weights must explicitly specify the units so that they produce a dimensionless quantity when multiplied by the residuals. 

#### Intermediate Simulation Outputs
- Rather than pass a `Val(true)` to the ODE RHS functions, to get intermediate calculated outputs, call the function `calc_md_Q` which dispatches on the parameter type to decide which model to use. This function returns a `NamedTuple` of the intermediate outputs, which can be used for plotting or analysis. The `Val(true)` argument is no longer accepted by the ODE RHS functions.

#### Other Breaking Changes
- The function `qrf_integrate` now returns a `NamedTuple` rather than a `Dict` with string keys.
- The fields `Arad` and `alpha` in `ParamObjRF` are removed, as they were unused.


### New Features

#### Modular Data Framework
- Introduced a robust hierarchy for experimental data:
    - `ExpFitData`: Container for multiple experimental datasets.
    - `AbstractExpDatum`: Base type for specific measurement types including `TfData`, `TvwSeriesData`, `TvwEndData`, and `EndTimeData`.
- Added `SolTrim` and `trim_sol` for precise alignment of model solutions with experimental time grids.

#### Extensibility Enhancements

### 🛠️ Refactors & Improvements
- Namespace Cleanup: Reduced namespace pollution by limiting re-exports.
- Error Handling: Improved validation in parameter constructors (e.g., checking units and dimensions).

## Previous Versions
A semi-formal changelog for previous versions can be found in the [GitHub releases](https://github.com/LyoHUB/LyoPronto.jl/releases).