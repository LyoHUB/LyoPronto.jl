# # Extending LyoPronto with New Models

# LyoPronto provides a flexible framework for lyophilization modeling. This document explains 
# how to implement a new model by walking through an extension of the Pikal model with 
# microwave heating as an example.

# The core pattern for adding a new model consists of five steps:

# 1. **Define a parameter container** (`ParamObj` subtype)
# 2. **Implement the ODE right-hand-side function**
# 3. **Provide an `ODEProblem` constructor** for your parameter type
# 4. **(Optional) Add `TransformVariables` transforms** for parameter fitting
# 5. **(Optional) Add plot recipes** for visualization

# Each step leverages Julia's multiple dispatch, so the fitting machinery in 
# [`gen_sol_pd`](@ref), [`obj_pd`](@ref), [`gen_nsol_pd`](@ref), and [`objn_pd`](@ref) 
# will work with your new model automatically—provided you follow the conventions below.

# ---

# # Imports

# The following imports are used throughout this document.

using LyoPronto
using TransformVariables
using OptimizationOptimJL
using LineSearches
using ADTypes: AutoForwardDiff
using Plots
using Accessors

# ---

# # 1. Define a parameter container

# Every model in LyoPronto has a parameter container type that extends the abstract type 
# `ParamObj`. This container holds all physical, geometric, and cycle parameters needed 
# by the ODE right-hand-side function.

# ### Conventions

# - The struct should be declared with `@concrete terse` (from `ConcreteStructs`). The 
#   `@concrete` macro ensures that all fields are concretely typed by inferring the most 
#   specific type from the values provided at construction time, rather than leaving them 
#   as `Any`. This is critical for performance: it avoids type instability throughout the 
#   ODE right-hand-side function and the fitting pipeline, where parameters are accessed 
#   repeatedly. The `terse` keyword suppresses the generation of default constructors, 
#   keeping the struct definition minimal.
# - Subtype `ParamObj` so the fitting machinery recognizes it.
# - Provide a constructor that accepts a tuple-of-tuples. This is the legacy format and 
#   is used by the fitting transforms.
# - Implement `Base.getindex` to return parameter groups as tuples (matching the tuple-of-tuples structure).
# - Implement `Base.size` to return the number of parameter groups.

# ### Example: `ParamObjPikalMW`

# The base Pikal model groups parameters into three tuples:

# ```julia
# params = (
#     (Rp, hf0, csolid, ρsolution),    # (1) Formulation parameters
#     (Kshf, Av, Ap),                  # (2) Heat transfer geometry
#     (pch, Tsh),                      # (3) Cycle controls
# )
# ```

# Our extended model adds a fourth group for microwave heating:

# ```julia
# params = (
#     (Rp, hf0, csolid, ρsolution),    # (1) Formulation parameters
#     (Kshf, Av, Ap),                  # (2) Heat transfer geometry
#     (pch, Tsh),                      # (3) Cycle controls
#     (P_mw, α_mw),                    # (4) Microwave heating (NEW)
# )
# ```

# The corresponding struct looks like:

# ```julia
# @concrete terse struct ParamObjPikalMW <: ParamObj
#     Rp
#     hf0
#     csolid
#     ρsolution
#     Kshf
#     Av
#     Ap
#     pch
#     Tsh
#     P_mw
#     α_mw
# end
# ```

# A tuple-of-tuples constructor and `Base.getindex` are also provided so that the fitting 
# machinery can access parameters by group index.

# ### For your model

# Design your parameter groups to match the physical meaning of your model. The number and 
# composition of groups is up to you—just be consistent between the struct, the constructor, 
# and `getindex`.

# ---

# # 2. Implement the ODE right-hand-side function

# The ODE function has the signature `func!(du, u, params, t)` where:

# - `du` is the output derivative vector (modified in-place).
# - `u` is the current state vector (unitless).
# - `params` is your `ParamObj` instance (or a tuple-of-tuples).
# - `t` is the current time (unitless, in hours).

# ### Unit conventions

# - **Time** `t` is unitless but represents hours. Dimensionalize inside the function: `tn = t * u"hr"`.
# - **State** `u` is unitless but has implicit units. Assign them internally: e.g. `Tf = u[2] * u"K"`.
# - **Derivatives** `du` should be stripped back to unitless values matching the implicit units: 
#   `du[2] = ustrip(u"K/hr", dTf)`.

# ### Example: `pikal_mw!`

# The extended Pikal model has two state variables: `[hf, Tf]` (remaining frozen layer 
# thickness, product temperature). It is implemented as a DAE (differential-algebraic 
# equation) with mass matrix `Diagonal([1.0, 0.0])`, where the algebraic constraint 
# enforces the pseudosteady-state energy balance.

# Following the Pikal model pattern, the physics is split into two functions:

# 1. A **helper function** (`calc_md_Q_mw`) that computes the physical quantities 
#    (mass flow, heat transfer terms) and returns them as a named tuple. This function 
#    can be reused independently for diagnostics or post-solution analysis.
# 2. The **RHS function** (`pikal_mw!`) that calls the helper, then writes derivatives 
#    to `du`.

# #### Helper: `calc_md_Q_mw`

# This function unpacks parameters, dimensionalizes state and time, computes all heat 
# and mass transfer terms, and returns a named tuple.

@inline function calc_md_Q_mw(u, po, t)
    (;Rp, hf0, csolid, ρsolution, Kshf, Av, Ap, pch, Tsh, P_mw, α_mw) = po
    
    td = t * u"hr"
    hf = u[1] * u"cm"
    Tf = u[2] * u"K"
    hd = hf0 - hf
    
    # Shelf heat transfer
    Q_shf = Kshf(pch(td)) * Av * (Tsh(td) - Tf) |> u"W"
    
    # Microwave heating (NEW term)
    Q_mw = P_mw * Ap * (1 - exp(-α_mw * hf)) |> u"W"
    
    # Sublimation mass transfer
    Tsub = Tf - Q_shf / k_ice / Ap * hf
    delta_p = calc_psub(Tsub) - pch(td)
    md = -Ap * delta_p / Rp(hd) |> u"g/hr"
    
    return (; md, Q_shf, Q_mw)
end

# #### RHS: `pikal_mw!`

# The RHS function calls the helper, then computes and writes the derivative and residuals.

function pikal_mw!(du, u, params, t)
    (;md, Q_shf, Q_mw) = calc_md_Q_mw(u, params, t)
    
    (; csolid, ρsolution, Ap) = params
    
    Q_sub = uconvert(u"W", md * ΔHsub)
    dhf_dt = min(0.0u"cm/hr", md / (ρsolution - csolid) / Ap |> u"cm/hr")
    
    du[1] = ustrip(u"cm/hr", dhf_dt)
    du[2] = ustrip(u"W", Q_sub + Q_shf + Q_mw)
end

# Because this is a DAE, wrap the RHS in an `ODEFunction` with a mass matrix:

const pikal_mw_f = ODEFunction{true, SciMLBase.AutoSpecialize}(
    pikal_mw!, mass_matrix=Diagonal([1.0, 0.0])
)

# ### For your model

# - Support both the named struct and tuple-of-tuples in your RHS function. This is required 
#   for compatibility with the fitting machinery, which constructs parameters via 
#   `setproperties` and may pass them as tuples.
# - If your model is a DAE (like the Pikal model), wrap the RHS in an `ODEFunction` with 
#   a `mass_matrix` argument, then pass the `ODEFunction` to `ODEProblem`.

# ---

# # 3. Provide an `ODEProblem` constructor

# The fitting functions [`gen_sol_pd`](@ref) and [`gen_nsol_pd`](@ref) call 
# `ODEProblem(param_obj)` to construct the problem. By defining a method for your 
# `ParamObj` subtype, you integrate with the entire fitting pipeline.

# ### Required helper methods

# Before defining `ODEProblem`, you need two helpers:

# #### `calc_u0(po::YourParamObj)`

# Returns the initial condition vector `u0` as a unitless `Vector{Float64}`.

function calc_u0(po::ParamObjPikalMW)
    return [ustrip(u"cm", po.hf0), ustrip(u"K", float(po.Tsh(0u"s")))]
end

# #### `get_tstops(po::YourParamObj)`

# Returns a sorted, unique vector of time points where cycle parameters change (ramp 
# transitions). These are used as `tstops` in the ODE solver for accuracy.

function get_tstops(po::ParamObjPikalMW)
    get_tstops((po.Tsh, po.pch))
end

# The generic `get_tstops` function on tuples calls `extract_ts` on each element, which 
# already handles `RampedVariable` and `LinearInterpolation` types.

# #### `ODEProblem(po::YourParamObj; u0=calc_u0(po), tspan=(0.0, 1000.0))`

function ODEProblem(po::ParamObjPikalMW; u0=calc_u0(po), tspan=(0.0, 1000.0))
    tstops = get_tstops(po)
    return ODEProblem{true, SciMLBase.FullSpecialize}(
        pikal_mw_f, u0, tspan, po;
        tstops=tstops, callback=end_drying_callback,
        initializealg=BrownFullBasicInit(), dt=0.1
    )
end

# Key points:
# - The first argument to `ODEProblem` is your RHS `ODEFunction` (or bare function for explicit ODEs).
# - `callback=end_drying_callback` terminates integration when drying is complete 
#   (frozen layer thickness approaches zero).
# - For DAEs, `initializealg=BrownFullBasicInit()` provides consistent initial conditions.
# - `tspan` should be generous—the callback will stop early.

# ---

# # 4. (Optional) Add `TransformVariables` transforms for fitting

# LyoPronto uses the `TransformVariables` package to map unconstrained optimization 
# parameters to physically meaningful ranges. Convenience functions already exist for 
# common parameters:

# - [`K_transform_basic`](@ref) — transforms for `Kshf` (shelf heat transfer coefficient)
# - [`Rp_transform_basic`](@ref) — transforms for `Rp` (product resistance)
# - [`KBB_transform_basic`](@ref) — transforms for wall/field heat transfer parameters
# - [`KBB_transform_bounded`](@ref) — bounded version using logistic transforms

# ### How transforms work

# A transform maps a flat `Vector{Float64}` to a `NamedTuple` of parameter values:

trans_K = K_transform_basic(5.0u"W/m^2/K")
# Maps a scalar to (; Kshf = ConstPhysProp(5.0u"W/m^2/K"))

# For multi-experiment fitting, combine transforms with `as`:

shared_trans = as((
    separate = as(Vector, trans_Rp, 3),  # 3 separate Rp sets
    shared = trans_K,                    # 1 shared Kv
))

# ### Adding transforms for new parameters

# To add a transform for a new parameter, compose `TVScale` and `TVExp` (or `TVLogistic`):

using TransformVariables: TVScale, TVExp
my_param_transform = as((; 
    myParam = TVScale(myParamGuess) ∘ TVExp() 
))

# The `TVExp()` ensures the parameter stays positive; `TVScale` sets the scale from a guess.

# ### How fitting uses transforms

# The fitting functions expect a tuple `tpf = (transform, param_objs, fitdats)`. 
# Internally they:

# 1. Call `transform(tr, fitlog)` to get a `NamedTuple` of parameters.
# 2. Call `setproperties(po, fitprm)` to merge fitted params into the base `ParamObj`.
# 3. Call `ODEProblem(new_po)` to construct and solve the ODE.
# 4. Compare the solution to `PrimaryDryFit` data.

# Because `setproperties` (from `ConstructionBase`) works on any struct, and `ODEProblem` 
# dispatches on your `ParamObj` subtype, **no additional code is needed** for fitting to work.

# ---

# # 5. (Optional) Add plot recipes

# LyoPronto uses `RecipesBase` to provide `plot` recipes for custom types. Two recipes 
# are particularly useful:

# #### `PrimaryDryFit` recipe

# Plots experimental data (product temperature vs. time) with markers.

# #### `RampedVariable` recipe

# Plots cycle parameter ramps (shelf temperature, chamber pressure, etc.).

# To add a recipe for your model's solution, define:

@recipe function f(::Type{Val{:pikal_mw_sol}}, sol)
    # Return plotting data from the ODE solution
    @series begin
        label := "Tf"
        sol.t, sol[:, 2]
    end
end

# ---

# # Summary checklist

# To add a new model to LyoPronto, implement:

# | Step | Required? | Purpose |
# |------|-----------|---------|
# | `ParamObj` subtype | Yes | Container for all model parameters |
# | Tuple-of-tuples constructor | Yes | Legacy format compatibility |
# | `Base.getindex` / `Base.size` | Yes | Parameter group access for fitting |
# | ODE RHS function | Yes | Model physics |
# | `ODEFunction` with mass matrix | Optional | For DAE models (e.g., Pikal-style) |
# | `calc_u0` method | Yes | Initial conditions |
# | `get_tstops` method | Yes | Solver accuracy at ramp transitions |
# | `ODEProblem` method | Yes | Problem construction for fitting |
# | `TransformVariables` transforms | Optional | Parameter fitting support |
# | Plot recipes | Optional | Visualization |

# Once these are in place, the existing functions 
# [`gen_sol_pd`](@ref), [`obj_pd`](@ref), [`gen_nsol_pd`](@ref), [`objn_pd`](@ref), 
# [`nls_pd`](@ref), and [`nls_pd!`](@ref) will work with your model automatically.
