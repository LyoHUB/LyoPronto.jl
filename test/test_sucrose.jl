vialsize = "10R"
rad_i, rad_o = get_vial_radii(vialsize)
Ap = π*rad_i^2  # cross-sectional area inside the vial
Av = π*rad_o^2 # vial bottom area

# Formulation parameters
c_solid = 0.05u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 0.8u"cm^2*Torr*hr/g"
A1 = 14u"cm*Torr*hr/g"
A2 = 1u"1/cm"
Rp = RpFormFit(R0, A1, A2)

# Cycle parameters
Vfill = 4u"mL" # ml
pch = RampedVariable(70u"mTorr")
T_shelf_final = (273.15 -25 )u"K"  # shelf temperature in K
T_shelf_0 = (273.15 - 45)u"K" # shelf temperature in K
ramp_rate = 1 *u"K/minute" # ramp rate deg/min
Tsh = RampedVariable([T_shelf_0, T_shelf_final], [ramp_rate], [])

KC = 3.58e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)

# Computed parameters based on above
hf0 = Vfill / Ap

params_bunch = [
    (Rp, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh)
]

# -------------------

tspan = (0.0, 100.0) # hours
u0 = [ustrip(u"cm", hf0), 233]
prob = ODEProblem(lyo_1d_dae_f, u0, tspan, tuple(params_bunch...))
sol = solve(prob, Rodas4P(autodiff=false), tstops=1/3, callback=LyoPronto.end_drying_callback)


# Results from Python
maxT = -32.1975
drytime = 45.8u"hr"

@test_broken sol.t[end]*u"hr" ≈ drytime rtol=0.01
@test sol.t[end]*u"hr" ≈ drytime rtol=0.1
@test sol[2,end]-273.15 ≈ maxT atol=0.1

