export lumped_cap_rf, lumped_cap_rf_julian!

# Constant properties
const rho_ice = 0.918u"g/cm^3" 
const rho_glass = 2.2u"g/cm^3"
const ΔH_sub = 678u"cal / g"
const e_0 = 8.854187e-12u"F/m" # permittivity of free space Coulomb^2/J/m
const σ = 5.670367e-8u"W/m^2/K^4" # Stefan-Boltzmann Constant
const k_sucrose = 0.139u"W/m/K"

function shapefac(Bi)
    charfunc(x) = -x*besselj1(x) + Bi*besselj0(x)
    λm = zeros(200)
    for i in eachindex(λm)
        if i == 1
            λm[i] = find_zero(charfunc, π/2)
        else
            λm[i] = find_zero(charfunc, λm[i-1]+π)
        end
    end
    Cm = @. 2/λm*besselj1(λm)/(besselj0(λm)^2 + besselj1(λm)^2)
    integ = sum(Cm .* besselj1.(λm) .* tanh.(λm))
end

const Bi_samp = 10.0 .^range(-2, 5, length=71)
const S_samp = shapefac.(Bi_samp)
const S_interp = linear_interpolation(Bi_samp, S_samp, extrapolation_bc=Line())


"""
    function lumped_cap_rf(du, u, p, t)

TODO: include all parameters
TODO: handle dimensionality
"""
function lumped_cap_rf(u, params, tn; energy_output=false)
    # R0, A1, A2, cSolid, ρ_solution = params[1]
    # K_shf_f = params[2]
    # h_f0, m_f0, A_p, cp_f, A_v, m_v, cp_v, A_rad= params[3]
    # shelf_ramp_func, pch, P_per_vial = params[4] # P_per_vial in W
    # f_RF, epp_f, epp_vw,  = params[5] # f_RF in 1/s
    # alpha, K_vwf, B_f, B_vw = params[6]
    Rp, h_f0, cSolid, ρ_solution = params[1]
    K_shf_f, A_v, A_p, = params[2]
    pch, Tsh, P_per_vial = params[3] 
    m_f0, cp_f, m_v, cp_v, A_rad = params[4]
    f_RF, epp_f, epp_vw = params[5]
    alpha, K_vwf, B_f, B_vw = params[6]
    
    t = tn*u"hr" # Dimensional time
    m_f = u[1]*u"g"
    T_f = u[2]*u"K"
    T_vw = u[3]*u"K"

    porosity = (ρ_solution - cSolid)/ρ_solution
    k_dry = k_sucrose*(1-porosity)

    rad = sqrt(A_p/π)
    V_vial = m_v / rho_glass

    K_shf = K_shf_f(pch(t))

    h_f = m_f/m_f0 * h_f0 # cm
    h_d = h_f0 - h_f
    A_sub = A_p + alpha*sqrt(abs(h_d))
    
    # R_p = R0 + A1*h_d/(1 + A2*h_d) #cm^2-s-Torr/g
    mflow = A_sub/Rp(h_d)*(calc_psub(T_f) - pch(t)) # g/s

    Bi = uconvert(NoUnits, K_vwf*rad/k_dry)
    Q_vwf = 2π*(K_vwf*rad*h_f + k_dry*(h_f0-h_f)*S_interp(Bi)) * (T_vw-T_f)

    Q_sub = mflow*ΔH_sub
    Q_shf = K_shf*A_v*(Tsh(t)-T_f)
    Q_shw = 0.9*σ*(Tsh(t)^4-T_vw^4)*A_rad # As if exchanging heat above vial?
    Qppp_RF_f = 2*pi*f_RF*e_0*epp_f*P_per_vial(t) * B_f # W / m^3
    Qppp_RF_vw = 2*pi*f_RF*e_0*epp_vw*P_per_vial(t) * B_vw # W / m^3
    Q_RF_f = Qppp_RF_f*A_p*h_f # W
    Q_RF_vw = Qppp_RF_vw*V_vial # W
    if Q_RF_f + Q_RF_vw > P_per_vial(t) && t == 0u"hr"
        @warn "Energy balance of EM terms not satisfied." uconvert(u"W", Q_RF_f) uconvert(u"W", Q_RF_vw) P_per_vial(t)
    end

    dm_f = min(0u"kg/s", -mflow/porosity)
    dT_f = (Q_RF_f -Q_sub + Q_vwf+ Q_shf) / (m_f * cp_f) - T_f*dm_f/m_f
    dT_vw = (Q_RF_vw + Q_shw - Q_vwf) / (m_v * cp_v)

    #     if Q_RF_f + Q_RF_vw > P_per_vial(t):
    #         print("Energy balance of EM terms not satisfied.")
    if energy_output
        # Q_vwf1 = 2π*(K_vwf*rad*h_f)*(T_vw-T_f)
        # Q_vwf2 = 2π*(k_dry*(h_f0-h_f))*(T_vw-T_f)
        return [ustrip(u"g/hr", dm_f), ustrip(u"K/hr", dT_f), ustrip(u"K/hr", dT_vw)], uconvert.(u"W", [Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw, ])
    else
        return [ustrip(u"g/hr", dm_f), ustrip(u"K/hr", dT_f), ustrip(u"K/hr", dT_vw)]
    end
end