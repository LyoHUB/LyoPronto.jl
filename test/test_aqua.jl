using Aqua
# TODO: after https://github.com/JuliaTesting/Aqua.jl/issues/387 or 
# https://github.com/JuliaLang/Pkg.jl/issues/4705 is resolved, can put persistent_tasks 
# check back in by deleting the kwarg
Aqua.test_all(LyoPronto; persistent_tasks=false)