using JET
# Can't use `mode = :basic` in CI because the Plots.jl recipe system triggers lots of false 
# positives for JET. Not a bad idea to periodically check that manually, though.
test_package(LyoPronto; mode = :typo, target_defined_modules=true)