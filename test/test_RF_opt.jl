using Optimization, OptimizationOptimJL
using LineSearches
optalg = LBFGS(linesearch=LineSearches.BackTracking())

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
# RF fit parameters (dummy values)
Bf = 2.0e7u"Ω/m^2"
Bvw = 0.9e7u"Ω/m^2"
Kvwf = 1.0e-3u"cal/s/K/cm^2"
# Controllable inputs
f_RF = 8u"GHz"
pch = RampedVariable(100u"mTorr")
Tsh = RampedVariable([233.15u"K", 283.15u"K"], 0.5u"K/minute",)
P_per_vial = RampedVariable(10u"W"/17 * 0.54) # actual power / vial

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
Tf = base_sol[2,:]*u"K"
Tvw = base_sol[3,:]*u"K"
t_end = t[end]
pdfit = PrimaryDryFit(t, Tf, Tvw, t_end)

tr = KBB_transform_basic(Kvwf*0.5, Bf*0.5, 0.5*Bvw)
pg = fill(1.0, 3)
sol = @inferred gen_sol_pd(pg, tr, po)
@test sol != base_sol
pass = (tr, po, pdfit)
# err = @inferred obj_pd(pg, pass)
err = @inferred obj_pd(pg, pass)
obj = OptimizationFunction(obj_pd, AutoForwardDiff(chunksize=3))
opt = solve(OptimizationProblem(obj, pg, pass), optalg;)
vals = transform(tr, opt.u)
@test vals.Kvwf ≈ Kvwf rtol=0.1
@test vals.Bf ≈ Bf rtol=0.5
@test vals.Bvw ≈ Bvw rtol=0.1


