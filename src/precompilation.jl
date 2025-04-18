@setup_workload begin
    @compile_workload begin

    # Conventional parameters
    # Geometry
    vialsize = "6R"
    rad_i, rad_o = get_vial_radii(vialsize)
    Ap = π*rad_i^2  # cross-sectional area inside the vial
    Av = π*rad_o^2 # vial bottom area
    Vfill = 3u"mL" # ml
    hf0 = Vfill / Ap
    # Formulation
    csolid = 0.05u"g/mL" # g solute / mL solution
    ρsolution = 1u"g/mL" # g/mL total solution density
    R0 = 0.8u"cm^2*Torr*hr/g"
    A1 = 14u"cm*Torr*hr/g"
    A2 = 1u"1/cm"
    Rp = RpFormFit(R0, A1, A2)
    # Cycle parameters
    pch = RampedVariable(70u"mTorr") # constant pressure
    T_shelf_0 = -40.0u"°C" # initial shelf temperature
    T_shelf_final = -10.0u"°C"  # final shelf temperature
    ramp_rate = 0.5 *u"K/minute" # ramp rate
    # Ramp for shelf temperature: convert to Kelvin because Celsius doesn't do math very well
    Tsh = RampedVariable(uconvert.(u"K", [T_shelf_0, T_shelf_final]), ramp_rate)
    # Heat transfer
    KC = 6.556e-5u"cal/s/K/cm^2"
    KP = 2.41e-3u"cal/s/K/cm^2/Torr"
    KD = 2.62u"1/Torr"
    Kshf = RpFormFit(KC, KP, KD)
    params_bunch = [
        (Rp, hf0, csolid, ρsolution),
        (Kshf, Av, Ap),
        (pch, Tsh)
    ]
    po = ParamObjPikal(params_bunch)

    # -------------------
    # Simulate conventional
    prob = ODEProblem(po)
    sol = solve(prob, Rosenbrock23())

    # ------------------- 
    # RF-specific parameters
    # Heat transfer
    cpf = cp_ice
    cpv = cp_gl
    mv = get_vial_mass(vialsize)
    mf0 = Vfill * ρsolution
    f_RF = 8u"GHz"
    # eppf = LyoPronto.ϵppf
    epp_w = epp_gl
    P_per_vial = RampedVariable(10u"W"/17 * 0.54)
    Bf = 2.0e7u"Ω/m^2"
    Bvw = 0.9e7u"Ω/m^2"
    Kvwf = 1.0e-3u"cal/s/K/cm^2"
    pch_rv = RampedVariable(100u"mTorr")
    Tsh_rv = RampedVariable([233.15u"K", 283.15u"K"], 0.5u"K/minute",)
    P_per_vial = RampedVariable(10u"W"/17 * 0.54) # actual power / vial

    params_bunch_rf = Vector{Any}([
        (Rp, hf0, csolid, ρsolution),
        (Kshf, Av, Ap),
        (pch_rv, Tsh_rv, P_per_vial),
        (mf0, cpf, mv, cpv),
        (f_RF, eppf, epp_w),
        (Kvwf, Bf, Bvw),
    ])
    po_rf = ParamObjRF(params_bunch_rf)

    # -------------------
    # Simulate RF
    prob = ODEProblem(po_rf)
    sol = solve(prob, Rosenbrock23())
    end
end