using Optimization, OptimizationOptimJL
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
T_shelf_0 = (273.15 -15 )u"K" # shelf temperature in K
T_shelf_final = (273.15 +10 )u"K"  # shelf temperature in K
ramp_rate = 0.5 *u"K/minute" # ramp rate deg/min
Tsh = RampedVariable([T_shelf_0, T_shelf_final], [ramp_rate], [])
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
pdfit = PrimaryDryFit(t, T, t_end)

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
    @test vals.Kshf(pch(0)) ≈ Kshf(pch(0)) rtol=0.3
    @test vals.Rp.R0 ≈ R0 rtol=0.1
    @test vals.Rp.A1 ≈ A1 rtol=0.1
    @test vals.Rp.A2 ≈ A2 rtol=0.1
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
    @test vals.Rp.A1 ≈ A1 rtol=0.1
    @test vals.Rp.A2 ≈ A2 rtol=0.5
end


@testset "Both Kv and Rp, least squares routine" begin
    tr = KRp_transform_basic(Kshf(pch(0))*0.75, R0*0.5, 2*A1, A2*0.5)
    pg = fill(0.0, 4)
    sol = @inferred gen_sol_pd(pg, tr, po)
    @test sol != base_sol
    pass = (tr, po, pdfit)
    nls = NonlinearFunction{true}(nls_pd, resid_prototype=zeros(num_errs(pdfit)))
    opt = solve(NonlinearLeastSquaresProblem(nls, pg, pass), GaussNewton(), reltol=1e-10, abstol=1e-10)
    vals = transform(tr, opt.u)
    @test vals.Kshf(pch(0)) ≈ Kshf(pch(0)) rtol=0.1
    @test vals.Rp.R0 ≈ R0 rtol=0.1
    @test vals.Rp.A1 ≈ A1 rtol=0.1
    @test vals.Rp.A2 ≈ A2 rtol=0.3
end