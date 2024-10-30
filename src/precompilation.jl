@setup_workload begin
    @compile_workload begin

    rad = 1.1u"cm"
    Ap = π*rad^2  # cross-sectional area inside the vial
    Av = π*(rad+1u"mm")^2 # vial bottom area

    # Formulation parameters
    c_solid = 0.05u"g/mL" # g solute / mL solution
    ρ_solution = 1u"g/mL" # g/mL total solution density
    R0 = 0.8u"cm^2*Torr*hr/g"
    A1 = 14u"cm*Torr*hr/g"
    A2 = 1u"1/cm"
    Rp = RpFormFit(R0, A1, A2)

    # Fill
    Vfill = 3u"mL" # ml
    hf0 = Vfill / Ap

    # Cycle parameters
    pch = RampedVariable(70u"mTorr") # constant pressure
    T_shelf_0 = -40.0u"°C" # initial shelf temperature
    T_shelf_final = -10.0u"°C"  # final shelf temperature
    ramp_rate = 0.5 *u"K/minute" # ramp rate
    # Ramp for shelf temperature: convert to Kelvin because Celsius doesn't do math very well
    Tsh = RampedVariable(uconvert.(u"K", [T_shelf_0, T_shelf_final]), ramp_rate)

    KC = 6.556e-5u"cal/s/K/cm^2"
    KP = 2.41e-3u"cal/s/K/cm^2/Torr"
    KD = 2.62u"1/Torr"
    Kshf = RpFormFit(KC, KP, KD)
    # Kshf = p-> 5.0u"W/m^2/K"

    params_bunch = [
        (Rp, hf0, c_solid, ρ_solution),
        (Kshf, Av, Ap),
        (pch, Tsh)
    ]

    # -------------------

    tspan = (0.0, 100.0) # hours

    u0 = [ustrip(u"cm", hf0), ustrip(u"K", Tsh(0u"minute"))]
    prob = ODEProblem(lyo_1d_dae_f, u0, tspan, tuple(params_bunch...))
    sol = solve(prob, Rosenbrock23(), callback=LyoPronto.end_drying_callback)



    end
end