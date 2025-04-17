vialsize = "6R"
rad_i, rad_o = get_vial_radii(vialsize)
Ap = π*rad_i^2  # cross-sectional area inside the vial
Av = π*rad_o^2 # vial bottom area
# Formulation parameters
c_solid = 0.06u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 0.8u"cm^2*Torr*hr/g"
A1 = 14u"cm*Torr*hr/g"
A2 = 1u"1/cm"
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
params_bunch = [
    (Rp, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh)
]
otherparams = (hf0, c_solid, ρ_solution, Av, Ap, pch, Tsh)
tspan = (0.0, 100.0) # hours
u0 = [ustrip(u"cm", hf0), 233]
KRp_prm = [10.0, 0.6, 12.0, 0.5]
po = ParamObjPikal(params_bunch)

# sol, new_params = gen_sol_conv_dim(KRp_prm, otherparams, u0, tspan)
# @test sol isa ODESolution
# @test size(new_params, 1) == 3
# @test new_params isa Tuple

# sol, new_params = gen_sol_conv_dim(KRp_prm, po, u0, tspan)
# @test sol isa ODESolution
# @test length(new_params) == 3
# @test new_params isa ParamObjPikal
# @test new_params.Rp isa RpFormFit
# @test new_params.Rp == RpFormFit(0.6*u"cm^2*Torr*hr/g", 12.0u"cm*Torr*hr/g", 0.5u"1/cm")
# @test ~(new_params.Kshf isa RpFormFit)
# @test new_params.Kshf(70u"mTorr") == 10u"W/m^2/K"
# @test new_params.Kshf(100u"mTorr") == 10u"W/m^2/K"