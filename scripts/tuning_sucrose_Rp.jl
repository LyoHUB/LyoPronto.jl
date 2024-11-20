
# Read in process data -------------- 

# procdat = CSV.File(datadir("exp_raw", "2024-06-21-16_MFD_AH.csv"), header=8)

# is_PD = procdat["Phase"] .== 4

# tproc = range(0, length=sum(is_PD), step=1/60)*u"hr"
# p_pir = procdat["VacPirani"][is_PD].*u"mTorr"
# Tsh_d = procdat["ShelfSetPT"][is_PD].*u"°C"
# T1 = procdat["TP1"][is_PD]u"°C"
# T2 = procdat["TP2"][is_PD]u"°C"
# T3 = procdat["TP4"][is_PD]u"°C"

# exptfplot(tproc, T1, T2, T3, labels=[L"T_{p1}" L"T_{p2}" L"T_{p3}"])
# plot!(tproc, Tsh_d, c=:black)

# T1trm = T1[tproc .< 7.5u"hr"]
# T2trm = T2[tproc .< 12.5u"hr"]
# T3trm = T3[tproc .< 7.5u"hr"]

# p_pir_sm = savitzky_golay(ustrip.(u"mTorr", p_pir), 91, 3, deriv=0).y *u"mTorr" # 91 a window width; 3 the polynomial order
# p_pir_der2 = savitzky_golay(ustrip.(u"mTorr", p_pir), 91, 3, deriv=2).y
# t_end = tproc[argmax(p_pir_der2[100:end])+99]
# plot(tproc, p_pir, label="data")
# plot!(tproc, p_pir_sm, label="smoothed data")
# tendplot!(t_end)


# Set up model --------------

# Vial geometry
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
sol = solve(prob, Rodas4P(), callback=LyoPronto.end_drying_callback)

plot(tproc, Tsh_d, c=:black, label=L"T_{sh}")
modconvtplot!(sol, label=L"$T_p$, model")


# ------------------ Match temperatures to get Rp, based on center vial
otherparams = (hf0, c_solid, ρ_solution, Av, Ap, pch, Tsh)
gen_sol_conv = KRp -> gen_sol_conv_dim(KRp, otherparams, u0, tspan)

obj_KRp1 = KRp -> obj_tT_conv(KRp, gen_sol_conv, (tproc[tproc.<7.5u"hr"], T1trm), t_end=t_end) 
obj_KRp2 = KRp -> obj_tT_conv(KRp, gen_sol_conv, (tproc[tproc.<12.5u"hr"], T2trm), t_end=t_end) 
obj_KRp3 = KRp -> obj_tT_conv(KRp, gen_sol_conv, (tproc[tproc.<7.5u"hr"], T3trm), t_end=t_end) 
obj_Rp1 = Rp -> obj_tT_conv([13, Rp...], gen_sol_conv, (tproc[tproc.<7.5u"hr"], T1trm), t_end=t_end) 
obj_Rp2 = Rp -> obj_tT_conv([13, Rp...], gen_sol_conv, (tproc[tproc.<12.5u"hr"], T2trm), t_end=t_end) 
obj_Rp3 = Rp -> obj_tT_conv([13, Rp...], gen_sol_conv, (tproc[tproc.<7.5u"hr"], T3trm), t_end=t_end) 

@time testconv = gen_sol_conv([22.5, 0.1, 10, 0.1])[1]
plot(testconv, idxs=2)

@time obj_KRp1([22.5, 0.1, 10, 0.1])
opt_KRp1 = optimize(obj_KRp1, [12, 0.1, 5, 0.1], NelderMead())
opt_KRp2 = optimize(obj_KRp2, [12, 0.1, 5, 0.1], NelderMead())
opt_KRp3 = optimize(obj_KRp3, [12, 0.1, 5, 0.1], NelderMead())
opt_Rp1 = optimize(obj_Rp1, [0.1, 10, 0.1], NelderMead())
opt_Rp2 = optimize(obj_Rp2, [0.1, 10, 0.1], NelderMead())
opt_Rp3 = optimize(obj_Rp3, [0.1, 10, 0.1], NelderMead())

KRp_1 = Optim.minimizer(opt_KRp1)
KRp_2 = Optim.minimizer(opt_KRp2)
KRp_3 = Optim.minimizer(opt_KRp3)
Rp_1 = Optim.minimizer(opt_Rp1)
Rp_2 = Optim.minimizer(opt_Rp2)
Rp_3 = Optim.minimizer(opt_Rp3)

prof1 = gen_sol_conv(KRp_1)[1]
prof2 = gen_sol_conv(KRp_2)[1]
prof3 = gen_sol_conv(KRp_3)[1]
profb1 = gen_sol_conv([13, Rp_1...])[1]
profb2 = gen_sol_conv([13, Rp_2...])[1]
profb3 = gen_sol_conv([13, Rp_3...])[1]


begin

exptfplot(tproc, T1, T2, T3)
modconvtplot!(prof1, prof2, prof3)
modconvtplot!(profb1, profb2, profb3, c=palette(:tab20c)[9:11])
plot!(ylim=(-42, -25), legend=:bottom, legend_columns=3)
end
