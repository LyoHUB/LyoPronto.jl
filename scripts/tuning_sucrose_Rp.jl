using DrWatson
@quickactivate :LyoPronto

using CSV
using SavitzkyGolay
using Optim
using LaTeXStrings
using Interpolations: linear_interpolation

# --------------- Set some plot defaults

default(:fontfamily, "Helvetica")
default(:palette, :tab20c)
default(:lw, 2)
resetfontsizes()
scalefontsizes(1.5)

# ------------
procdat = CSV.File(datadir("exp_raw", "2024-06-21-16_MFD_AH.csv"), header=8)

# procdat["Phase"]
is_PD = procdat["Phase"] .== 4

tproc = range(0, length=sum(is_PD), step=1/60)*u"hr"
p_pir = procdat["VacPirani"][is_PD].*u"mTorr"
Tsh_d = procdat["ShelfSetPT"][is_PD].*u"°C"
plot(p_pir)

p_pir_sm = savitzky_golay(ustrip.(u"mTorr", p_pir), 91, 3, deriv=0).y *u"mTorr"
p_pir_der2 = savitzky_golay(ustrip.(u"mTorr", p_pir), 91, 3, deriv=2).y
plot(tproc[100:end], p_pir_der2[100:end])
t_end = tproc[argmax(p_pir_der2[100:end])+99]
plot(tproc, p_pir)
plot!(tproc, p_pir_sm)
vline!([t_end])

T1 = procdat["TP1"][is_PD]u"°C"
T2 = procdat["TP2"][is_PD]u"°C"
T3 = procdat["TP4"][is_PD]u"°C"

plot(tproc, T1)
plot!(tproc, T2)
plot!(tproc, T3)

T1trm = T1[tproc .< 7.5u"hr"]
T2trm = T2[tproc .< 12.5u"hr"]
T3trm = T3[tproc .< 7.5u"hr"]

# --------------

ΔHsub = 678u"cal/g"

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

# Cycle parameters
Vfill = 3u"mL" # ml
pch = RampedVariable(70u"mTorr")
T_shelf_final = (273.15 -10 )u"K"  # shelf temperature in K
T_shelf_0 = (273.15 - 40)u"K" # shelf temperature in K
ramp_rate = 0.5 *u"K/minute" # ramp rate deg/min
Tsh = RampedVariable([T_shelf_0, T_shelf_final], [ramp_rate], [])

KC = 6.556e-5u"cal/s/K/cm^2"
KP = 2.41e-3u"cal/s/K/cm^2/Torr"
KD = 2.62u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)
# Kshf = p-> 5.0u"W/m^2/K"


# Computed parameters based on above
hf0 = Vfill / Ap
params_bunch = [
    (Rp, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh)
]

# -------------------

tspan = (0.0, 100.0) # hours
end_cond(u, t, integ) = u[1] 
end_drying = ContinuousCallback(end_cond, terminate!)

u0 = [ustrip(u"cm", hf0), 233]
prob = ODEProblem(lyo_1d_dae_f, u0, tspan, tuple(params_bunch...))
sol = solve(prob, Rodas4P(), callback=end_drying)

plot(sol, idxs=1)
plot(sol, idxs=2)

# -----------------

# T_sh_interp = linear_interpolation(tproc, uconvert.(u"K", Tsh_d))
# plot(tproc, Tsh_d)
# # plot!(t_idc, T_sh_interp.(t_idc))
# plot!(xlim=(0, 1))

# dT = T_sh_interp.(t_idc) - uconvert.(u"K", T_idc  )
# integ_dT = sum((dT[1:end-1] .+ dT[2:end])./2 .* (t_idc[2:end] .- t_idc[1:end-1]))

# Kv_est = ΔHsub*Vfill*(ρ_solution - c_solid) / Av / integ_dT

# uconvert(u"W/m^2/K", Kv_est)
# uconvert(u"cal/s/cm^2/K", Kv_est)

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
plot()
plot!(tproc, T1, c=1)
plot!(tproc, T2, c=2)
plot!(tproc, T3, c=3)
plot!(prof1.t.*u"hr", prof1.(prof1.t, idxs=2)*u"K", c=5)
plot!(prof2.t.*u"hr", prof2.(prof2.t, idxs=2)*u"K", c=6)
plot!(prof3.t.*u"hr", prof3.(prof3.t, idxs=2)*u"K", c=7)
# plot!(profb1.t, profb1.(profb1.t, idxs=2)*u"K", c=9)
# plot!(profb2.t, profb2.(profb2.t, idxs=2)*u"K", c=10)
# plot!(profb3.t, profb3.(profb3.t, idxs=2)*u"K", c=11)
plot!(ylim=(-42, -25), legend=:bottom, legend_columns=3)
end
