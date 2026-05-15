using TransformVariables
using OptimizationOptimJL
using LineSearches
using NonlinearSolve
optalg = LBFGS(linesearch=LineSearches.BackTracking())

vialsize = "6R"
rad_i, rad_o = get_vial_radii(vialsize)
Ap = π*rad_i^2  # cross-sectional area inside the vial
Av = π*rad_o^2 # vial bottom area
# Formulation parameters
csolid = 0.06u"g/mL" # g solute / mL solution
ρsolution = 1u"g/mL" # g/mL total solution density
R0 = 0.8u"cm^2*Torr*hr/g"
A1 = 14.0u"cm*Torr*hr/g"
A2 = 1.0u"1/cm"
Rp = RpFormFit(R0, A1, A2)
# Cycle parameters
Vfill = 3u"mL" # ml
pch = RampedVariable(70u"mTorr")
Tsh = RampedVariable([-15u"°C", 10u"°C"].|>u"K", 0.5u"K/minute")
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)
# Computed parameters based on above
hf0 = Vfill / Ap
po = ParamObjPikal((
    (Rp, hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch, Tsh)
))
base_sol = solve(ODEProblem(po), LyoPronto.odealg_chunk2)

t = base_sol.t*u"hr"
T = base_sol[2,:]*u"K"
t_end = t[end]
pdfit = PrimaryDryFit(t, T; t_end)

@testset "Both Kv and Rp, optimization routine" begin
    tr = KRp_transform_basic(Kshf(pch(0))*0.75, R0*0.5, 2*A1, A2*0.5)
    pg = fill(0.0, 4)
    sol = @inferred gen_sol_pd(pg, tr, po)
    @test sol != base_sol
    pass = (tr, po, pdfit)
    # err = @inferred obj_pd(pg, pass)
    err = @inferred obj_pd(pg, pass)
    obj = OptimizationFunction(obj_pd, AutoForwardDiff(chunksize=4))
    opt = solve(OptimizationProblem(obj, pg, pass), optalg)
    vals = transform(tr, opt.u)
    @test all(opt.u .!= 0)
    @test vals.Kshf(pch(0)) ≈ Kshf(pch(0)) rtol=0.3
    @test vals.Rp.R0 ≈ R0 rtol=0.1
    @test vals.Rp.A1 ≈ A1 rtol=0.3
    @test vals.Rp.A2 ≈ A2 rtol=0.5
end

@testset "Only Rp" begin
    tr = Rp_transform_basic(R0*0.5, 2*A1, A2*0.5)
    pg = fill(0.0, 3)
    sol = @inferred gen_sol_pd(pg, tr, po)
    @test sol != base_sol
    pass = (tr, po, pdfit)
    # err = @inferred obj_pd(pg, pass)
    err = @inferred obj_pd(pg, pass)
    obj = OptimizationFunction(obj_pd, AutoForwardDiff(chunksize=3))
    opt = solve(OptimizationProblem(obj, pg, pass), optalg)
    vals = transform(tr, opt.u)
    @test vals.Rp.R0 ≈ R0 rtol=0.1
    @test vals.Rp.A1 ≈ A1 rtol=0.2
    @test vals.Rp.A2 ≈ A2 rtol=0.5

    # Check that the badprms path outputs NaNs as expected
    badprms = x->true
    @test isnan(obj_pd(pg, pass; badprms))

    # Check that it also works if given a specific condition
    badprms = x->x.Kshf(pch(0)) < 100u"W/m^2/K" 
    @test isnan(obj_pd(pg, pass; badprms))
end


@testset "Both Kv and Rp, least squares routine" begin
    tr = KRp_transform_basic(Kshf(pch(0))*0.75, R0*0.5, 2*A1, A2*0.5)
    pg = fill(0.0, 4)
    sol = @inferred gen_sol_pd(pg, tr, po)
    @test sol != base_sol
    pass = (tr, po, pdfit)
    nls = NonlinearFunction(pdfit)
    opt = solve(NonlinearLeastSquaresProblem(nls, pg, pass), GaussNewton(), reltol=1e-10, abstol=1e-10)
    vals = transform(tr, opt.u)
    @test vals.Kshf(pch(0)) ≈ Kshf(pch(0)) rtol=0.1
    @test vals.Rp.R0 ≈ R0 rtol=0.1
    @test vals.Rp.A1 ≈ A1 rtol=0.1
    @test vals.Rp.A2 ≈ A2 rtol=0.3
end

po2 = @set po.Rp = RpFormFit(2u"cm^2*Torr*hr/g", 5u"cm*Torr*hr/g", 1.5u"cm^-1")
po3 = @set po.Rp = RpFormFit(0.5u"cm^2*Torr*hr/g", 20u"cm*Torr*hr/g", 0.0u"cm^-1")

pos = [po, po2, po3]
pdfits = map(pos) do poi
    base_sol = solve(ODEProblem(poi), LyoPronto.odealg_chunk2)
    t = base_sol.t*u"hr"
    T = base_sol[2,:]*u"K"
    t_end = t[end]
    pdfit = PrimaryDryFit(t, T; t_end)
end


@testset "Fit with shared Kv, distinct Rp" begin
    big_trans = as((; 
        shared = K_transform_basic(Kshf(pch(0))*0.75),
        separate = as(Vector, Rp_transform_basic(R0*0.75, A1*2, A2*0.5), 3)
    ))



    pg = fill(0.0, TransformVariables.dimension(big_trans))
    @test_broken @inferred gen_nsol_pd(pg, big_trans, pos)
    sols = gen_nsol_pd(pg, big_trans, pos)
    @testset "Different solution" begin
        for sol in sols
            @test sol != base_sol
        end
    end
    pass = (big_trans, pos, pdfits)
    # err = @inferred objn_pd(pg, pass)

    # This specific test can be deleted if it becomes trouble, probably
    exact = log.([1/0.75, 1/0.75, 0.5, 2, 2/0.8/0.75, 5/14/2, 1.5*2, 0.5/0.8/0.75, 20/14/2, 1e-20])
    @test objn_pd(exact, pass) ≈ 0 atol=1e-4  # Transformation should give zero objective at original values

    obj = OptimizationFunction(objn_pd, AutoForwardDiff(chunksize=5)) # Length of 10: divide it nicely
    opt = solve(OptimizationProblem(obj, pg, pass), optalg, f_abstol=1e-2)
    vals = transform(big_trans, opt.u)
    @test vals.shared.Kshf(pch(0)) ≈ Kshf(pch(0)) rtol=0.1
    for (poi, sepi) in zip(pos, vals.separate)
        @test sepi.Rp.R0 ≈ poi.Rp.R0 rtol=0.1
        @test sepi.Rp.A1 ≈ poi.Rp.A1 rtol=0.3
        @test sepi.Rp.A2 ≈ poi.Rp.A2 atol=0.6u"cm^-1"
    end

    # Check that the badprms path outputs NaNs as expected
    badprms = x->true
    @test all(isnan.(objn_pd(pg, pass; badprms)))

    badprms = x->x.Kshf(pch(0)) < 100u"W/m^2/K" 
    @test all(isnan.(objn_pd(pg, pass; badprms)))
end
