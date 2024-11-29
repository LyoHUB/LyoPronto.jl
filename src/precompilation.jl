@setup_workload begin
    @compile_workload begin

    # Conventional parameters
    # Geometry
    rad = 1.1u"cm"
    Vfill = 3u"mL" # ml
    Ap = π*rad^2  # cross-sectional area inside the vial
    Av = π*(rad+1u"mm")^2 # vial bottom area
    hf0 = Vfill / Ap
    # Formulation
    c_solid = 0.05u"g/mL" # g solute / mL solution
    ρ_solution = 1u"g/mL" # g/mL total solution density
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
        (Rp, hf0, c_solid, ρ_solution),
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
    cp_f = 2.050u"J/g/K"
    cp_v = 0.840u"J/g/K"
    m_v = get_vial_mass(vialsize)
    m_f0 = Vfill * ρ_solution
    f_RF = 8u"GHz"
    epp_f = LyoPronto.ϵpp_f
    epp_w = 2.4e-2
    P_per_vial = RampedVariable(10u"W"/17 * 0.54)
    B_f = 2.0e7u"Ω/m^2"
    B_vw = 0.9e7u"Ω/m^2"
    K_vwf = 1.0e-3u"cal/s/K/cm^2"
    pch_rv = RampedVariable(100u"mTorr")
    Tsh_rv = RampedVariable([233.15u"K", 283.15u"K"], 0.5u"K/minute",)
    P_per_vial = RampedVariable(10u"W"/17 * 0.54) # actual power / vial

    params_bunch_rf = Vector{Any}([
        (Rp, h_f0, c_solid, ρ_solution),
        (K_shf_f, A_v, A_p),
        (pch_rv, Tsh_rv, P_per_vial),
        (m_f0, cp_f, m_v, cp_v),
        (f_RF, epp_f, epp_w),
        (K_vwf, B_f, B_vw),
    ])
    po_rf = ParamObjRF(params_bunch_rf)

    # -------------------
    # Simulate RF
    prob = ODEProblem(po_rf)
    sol = solve(prob, Rosenbrock23(autodiff=false))
    end
end