using Optimization, OptimizationOptimJL
using LineSearches
optalg = LBFGS(linesearch=LineSearches.BackTracking())

vialsize = "6R"
rad_i, rad_o = get_vial_radii(vialsize)
A_p = π*rad_i^2  # cross-sectional area inside the vial
A_v = π*rad_o^2 # vial bottom area
m_v = get_vial_mass(vialsize)
# Formulation and fill
c_solid = 0.05u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16.0u"cm*hr*Torr/g"
A2 = 0.0u"1/cm"
Rp = RpFormFit(R0, A1, A2)
Vfill = 5u"mL"
# Heat transfer
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
K_shf_f = RpFormFit(KC, KP, KD)
# Geometry
h_f0 = Vfill/A_p
m_f0 = Vfill * ρ_solution
# RF fit parameters (dummy values)
B_f = 2.0e7u"Ω/m^2"
B_vw = 0.9e7u"Ω/m^2"
K_vwf = 1.0e-3u"cal/s/K/cm^2"
# Controllable inputs
f_RF = 8u"GHz"
pch = RampedVariable(100u"mTorr")
Tsh = RampedVariable([233.15u"K", 283.15u"K"], 0.5u"K/minute",)
P_per_vial = RampedVariable(10u"W"/17 * 0.54) # actual power / vial

po = ParamObjRF((
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf_f, A_v, A_p),
    (pch, Tsh, P_per_vial),
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl),
    (f_RF, LyoPronto.epp_f, LyoPronto.epp_gl),
    (K_vwf, B_f, B_vw),
))

base_sol = solve(ODEProblem(po), Rodas3())

t = base_sol.t*u"hr"
Tf = base_sol[2,:]*u"K"
Tvw = base_sol[3,:]*u"K"
t_end = t[end]
pdfit = PrimaryDryFit(t, T, Tvw, t_end)

tr = KBB_transform_basic(K_vwf*0.75, B_f*0.5, 2*B_vw)
pg = fill(0.0, 3)
sol = gen_sol_pd(pg, tr, po)
@test sol != base_sol
pass = (tr, po, pdfit)
# err = @inferred obj_pd(pg, pass)
err = obj_pd(pg, pass)
obj = OptimizationFunction(obj_pd, AutoForwardDiff())
opt = solve(OptimizationProblem(obj, pg, pass), optalg)
vals = transform(tr, opt.u)
@test vals.K_vwf ≈ K_vwf rtol=0.1
@test vals.B_f ≈ B_f rtol=0.1
@test vals.B_vw ≈ B_vw rtol=0.1

