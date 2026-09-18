using TransformVariables
using OptimizationOptimJL
using NonlinearSolve
using LineSearches
using Unitful
optalg = Optim.BFGS(linesearch=LineSearches.BackTracking())

vialsize = "6R"
rad_i, rad_o = get_vial_radii(vialsize)
Ap = π*rad_i^2  # cross-sectional area inside the vial
Av = π*rad_o^2 # vial bottom area
mv = get_vial_mass(vialsize)
# Formulation and fill
csolid = 0.05u"g/mL" # g solute / mL solution
ρsolution = 1u"g/mL" # g/mL total solution density
R0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16.0u"cm*hr*Torr/g"
A2 = 0.0u"1/cm"
Rp = RpFormFit(R0, A1, A2)
Vfill = 5u"mL"
# Heat transfer
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf_f = RpFormFit(KC, KP, KD)
# Geometry
hf0 = Vfill/Ap
mf0 = Vfill * ρsolution
# RF fit parameters (base to which we will fit)
Bf = 5.0e8u"Ω/m^2"
Bvw = 9.0e6u"Ω/m^2"
Kvwf = 20.0u"W/K/m^2"
# Controllable inputs
f_RF = 8u"GHz"
pch = RampedVariable(100u"mTorr")
Tsh = RampedVariable([233.15u"K", 283.15u"K"], 0.5u"K/minute",)
P_per_vial = RampedVariable(0.5u"W") 

po = ParamObjRF((
    (Rp, hf0, csolid, ρsolution),
    (Kshf_f, Av, Ap),
    (pch, Tsh, P_per_vial),
    (mf0, LyoPronto.cp_ice, mv, LyoPronto.cp_gl),
    (f_RF, LyoPronto.eppf, LyoPronto.epp_gl),
    (Kvwf, Bf, Bvw),
))

base_sol = solve(ODEProblem(po), LyoPronto.odealg_chunk3)

t = base_sol.t*u"hr"
keep = findall(diff(base_sol.t) .> .01)
t = base_sol.t[keep]*u"hr"
Tf = base_sol[2,keep]*u"K"
Tvw = base_sol[3,keep]*u"K"
t_end = t[end]
pdfit = PrimaryDryFit(t, Tf, Tvw, t_end)

tr = KBB_transform_basic(Kvwf*0.5, Bf*0.5, 0.5*Bvw)
pg = fill(1.0, 3)
sol = @inferred gen_sol_pd(pg, tr, po)
@test sol != base_sol
pass = (tr, po, pdfit)

@testset "qrf_integrate" begin
    qinteg = qrf_integrate(base_sol, po)

    @test qinteg isa Dict
    @test haskey(qinteg, "Qsub")
    @test haskey(qinteg, "Qshf")
    @test haskey(qinteg, "Qvwf")
    @test haskey(qinteg, "QRFf")
    @test haskey(qinteg, "QRFvw")
    @test haskey(qinteg, "Qshw")

    for v in values(qinteg)
        @test first(v) isa Unitful.Energy
    end

    # @test qinteg["Qsub"] > 0u"W*hr"

    # Energy conservation check on the product:
    #   d(mf*cpf*Tf)/dt = Q_shf + Q_vwf + Q_RF_f - Q_sub
    # Integrating:  ∫(Q_shf+Q_vwf+Q_RF_f) dt = ∫Q_sub dt + Δ(mf*cpf*Tf)
    m_f = base_sol[1, :] * u"g"
    T_f = base_sol[2, :] * u"K"
    cpf = po.cpf

    Δinternal = cpf * (m_f[end] * T_f[end] - m_f[begin] * T_f[begin]) |> u"W*hr"
    energy_in = qinteg["Qshf"] + qinteg["Qvwf"] + qinteg["QRFf"]
    energy_out = qinteg["Qsub"] + Δinternal

    @test isapprox(energy_in, energy_out; rtol=1e-2)

    # Energy conservation check on the vial wall:
    #   dTvw/dt = (Q_shw - Q_vwf + Q_RF_vw) / (mv*cpv)
    #   => mv*cpv*dTvw/dt = Q_shw - Q_vwf + Q_RF_vw
    # Integrating:  ∫(Q_shw - Q_vwf + Q_RF_vw) dt = Δ(mv*cpv*Tvw)
    T_vw = base_sol[3, :] * u"K"
    mv = po.mv
    cpv = po.cpv

    Δinternal_vw = mv * cpv * (T_vw[end] - T_vw[begin]) |> u"W*hr"
    energy_in_vw = qinteg["Qshw"] + qinteg["QRFvw"]
    energy_out_vw = qinteg["Qvwf"] + Δinternal_vw

    @test isapprox(energy_in_vw, energy_out_vw; rtol=1e-2)
end


@testset "Optimization" begin
    err = @inferred obj_pd(pg, pass)
    obj = OptimizationFunction(obj_pd, AutoForwardDiff(chunksize=3))
    opt = solve(OptimizationProblem(obj, pg, pass), optalg;)
    @test SciMLBase.successful_retcode(opt)
    vals = transform(tr, opt.u)
    @test vals.Kvwf ≈ Kvwf rtol=0.1
    @test vals.Bf ≈ Bf rtol=0.1
    @test vals.Bvw ≈ Bvw rtol=0.1
end

@testset "Least squares" begin
    lsq = NonlinearFunction(pdfit)
    opt = @inferred solve(NonlinearLeastSquaresProblem(lsq, pg, pass), LevenbergMarquardt())
    @test SciMLBase.successful_retcode(opt)
    vals = transform(tr, opt.u)
    @test vals.Kvwf ≈ Kvwf rtol=0.2
    @test vals.Bf ≈ Bf rtol=0.2
    @test vals.Bvw ≈ Bvw rtol=0.2
end


