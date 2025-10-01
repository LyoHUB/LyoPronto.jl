module LyoPronto

using Reexport
@reexport using OrdinaryDiffEqRosenbrock
import OrdinaryDiffEqRosenbrock: ODEProblem
@reexport using OrdinaryDiffEqNonlinearSolve
@reexport using DiffEqCallbacks
@reexport using Unitful: @u_str, ustrip
import Unitful
using TransformVariables
using RecipesBase
using ColorTypes: RGB
using CSV
using PrecompileTools
using SpecialFunctions: besselj0, besselj1
using DataInterpolations
using SavitzkyGolay
using Roots
using Accessors
using ConcreteStructs
using ADTypes: AutoForwardDiff
@reexport import ConstructionBase: setproperties
using DocStringExtensions
using LinearAlgebra: Diagonal

abstract type ParamObj end

const odealg_chunk1 = Rodas4(autodiff=AutoForwardDiff(chunksize=1))
const odealg_chunk2 = Rodas4(autodiff=AutoForwardDiff(chunksize=2))
const odealg_chunk3 = Rodas4(autodiff=AutoForwardDiff(chunksize=3))

include("structs.jl")
include("rf_lumcap_model.jl")
include("pikal_model.jl")
include("paramfits.jl")
include("recipes.jl")
include("cycle_time.jl")
include("get_vial_dims.jl")
include("physical_properties.jl")

# Exports, all in one place
# convenience structs
export RpFormFit, RampedVariable, ConstPhysProp, PrimaryDryFit
# simulation helpers
export end_drying_callback
export calc_u0, get_tstops
# conventional lyo
export lyo_1d_dae_f
export ParamObjPikal
export RpEstimator, calc_hRp_T
# RF lyo
export lumped_cap_rf!
export ParamObjRF
# parameter fitting tools
export gen_sol_pd, obj_pd, gen_nsol_pd, objn_pd
export KRp_transform_basic, K_transform_basic, Rp_transform_basic, KBB_transform_basic
export obj_expT, err_expT, err_expT!, num_errs, nls_pd, nls_pd!
export ConstWrapTV
# plotting tools (mostly already done by macros)
export qrf_integrate
# End of primary drying
export identify_pd_end
# Vial dimensions
export get_vial_radii, get_vial_mass, get_vial_shape, make_outlines


# include("precompilation.jl")

# If on 1.11 or later, mark stuff as public that isn't exported
@static VERSION >= v"1.11" && include("public.jl")

end