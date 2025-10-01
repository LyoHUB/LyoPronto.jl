module LyoPronto

using Reexport
@reexport using OrdinaryDiffEqRosenbrock
import OrdinaryDiffEqRosenbrock: ODEProblem
@reexport using OrdinaryDiffEqNonlinearSolve
@reexport using DiffEqCallbacks
@reexport using Unitful
@reexport using TransformVariables
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
# using .Dielectric

# include("precompilation.jl")

# If on 1.11 or later, mark stuff as public that isn't exported
VERSION >= v"1.11" && include("public.jl")

end