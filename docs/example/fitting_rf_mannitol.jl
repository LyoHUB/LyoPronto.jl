# # Imports
# LyoPronto is this package. It reexports several other packages, so after
# `using LyoPronto`, you have effectively also done `using Unitful` and a few others.
using LyoPronto 

## CSV and TypedTables load and store experimental data
using TypedTables, CSV

## Plots is a frontend to several plotting packages, and defaults to GR
using Plots
# Some plot defaults, which will apply to subsequent plots in a script/session
default(:fontfamily, "Computer Modern")
default(:framestyle, :box)
default(:lw, 2)
default(:markersize, 4)
default(:markerstrokewidth, 0.5)
default(:unitformat, :square)
resetfontsizes()
scalefontsizes(1.2)

## Latexify helps with LaTeX formatting in plots
using Latexify
## 
set_default(labelformat=:square) # Latexify, not Plots



# ----------------------------------------
# Load processed data into memory
reloaded = load(datadir("exp_pro", "M1_processed.jld2"))
@unpack thm_pd, lyo_pd, fitdat = reloaded

begin
# --------------------------
# - Set up parameters
# Vial parameters
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
Bf = 2.0e7u"Ω/m^2"
Bvw = 0.9e7u"Ω/m^2"
Kvwf = 10.0u"W/K/m^2"
# Controllable inputs
f_RF = 8u"GHz"
pch = RampedVariable(100u"mTorr")
Tsh = RampedVariable([233.15u"K", 283.15u"K"], 0.5u"K/minute",)
P_per_vial = RampedVariable(10u"W"/17 * 0.54) # actual power / vial

params_base = ParamObjRF((
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf_f, A_v, A_p),
    (pch, Tsh, P_per_vial),
    (m_f0, LyoPronto.cp_ice, m_v, LyoPronto.cp_gl),
    (f_RF, LyoPronto.eppf, LyoPronto.epp_gl),
    (Kvwf, Bf, Bvw),
))
end

# --------------- Perform the optimization
trans_KBB = KBB_transform_bounded(Kvwf, Bf, Bvw)
nls_M1 = NonlinearFunction{true}(nls_pd!, resid_prototype=zeros(num_errs(fitdat)))
p0 = [3.0, 3.0, 0.3]
tsol = gen_sol_pd(p0, trans_KBB, params_base)
modrftplot(tsol)
plot!(fitdat, nmarks=40)


# @time objf_KBB(p0, (trans_KBB, params_base, fitdat))
# opt1 = solve(OptimizationProblem(objf_KBB, p0, (trans_KBB, params_base, fitdat)), optalg, maxiters=100, show_trace=true)

opt1 = solve(NonlinearLeastSquaresProblem(nls_M1, p0, (trans_KBB, params_base, fitdat)), LevenbergMarquardt())
prof_RF = gen_sol_pd(opt1.u, trans_KBB, params_base)
transform(trans_KBB, opt1.u)

# Save fit results to a file
save_fitresults(opt1, "M1")
default(:linewidth, 3)

resetfontsizes()
scalefontsizes(1.2)
begin
plT = blankplot_hrC(margin_left=200Plots.px)
@df thm_pd exptfplot!(:t, :T4, :T1, nmarks=40) 
@df thm_pd exptvwplot!(:t, :T3, nmarks=40, label=L"$T_\mathrm{vw1}$, exp.")
modrftplot!(prof_RF, markeralpha=0, trimend=1)
plot!(Tsh, c=:black, label=L"T_\mathrm{sh}")
tendplot!(fitdat.t_end, label="", ls=:dash)
# tendplot!(thm_pd.t[argmax(thm_pd.T4)], ls=:dash, label="RF off")
plot!(legend=:topleft, ylim=(-40, 50), xlim=(0,13))
annotate!(8, -32, Plots.text("end of drying,\nRF off", 12, "Computer Modern"))
plot!([fitdat.t_end-1.5u"hr", fitdat.t_end-0.2u"hr"], [-30, -30], arrow=:arrow, c=:black, linewidth=1, label="")
plot!(size=(480,400), left_margin=20Plots.px)
# end
# # savefig(plotsdir("M1_compT.svg"))
# # savefig(plotsdir("M1_compT.pdf"))

# begin
plq = qplotrf(prof_RF, ordering=3:-1:1)
plot!(size=(400,300), ylim=(0,0.5), widen=false, left_margin=20Plots.px, bottom_margin=20Plots.px)
# savefig(plotsdir("M1_energy_budget.svg"))
# savefig(plotsdir("M1_energy_budget.pdf"))

plot(plT, plq, layout=@layout([a  b{0.4w}]), size=(800,400))
end
savefig(plotsdir("M1_T_q_combine.svg"))
savefig(plotsdir("M1_T_q_combine.pdf"))

begin
plT = blankplot_hrC(margin_left=200Plots.px)
@df thm_pd exptfplot!(:t, :T4, nmarks=40, label="exp, f") 
@df thm_pd exptvwplot!(:t, :T3, nmarks=40, msw=4, label="exp, vw")
modrftplot!(prof_RF, trimend=1, linealpha=0.8, lw=4, label=["LC: f, vw "  ""])
plot!(Tsh, c=:black, label="shelf")
tendplot!(fitdat.t_end, label="", ls=:dash)
# tendplot!(thm_pd.t[argmax(thm_pd.T4)], ls=:dash, label="RF off")
plot!(legend=:topleft, ylim=(-40, 50), xlim=(0,13))
annotate!(7.5, -28, Plots.text("end of drying,\nRF off", 12, "Computer Modern"))
plot!([fitdat.t_end-1.5u"hr", fitdat.t_end-0.2u"hr"], [-30, -30], arrow=:arrow, c=:black, linewidth=1, label="")
plot!(size=(360,250), left_margin=20Plots.px, legend_columns=2, legend=:topleft, legendfontsize=10)
plq = qplotrf(prof_RF, tot_lab="")
plot!(size=(400,300), ylim=(0,0.5), widen=false, left_margin=20Plots.px, bottom_margin=20Plots.px)
plot!(plq, legendfontsize=11, )
plot(plT, plq, layout=@layout([a  b{0.4w}]), size=(600,250))
end

savefig(plotsdir("M1_T_q_combine_poster.svg"))

begin
blankplot_hrC(size=(400,300), legend=false)
plot!(Tsh, c=:gray, tmax=13.9u"hr", label="")
@df thm_pd exptfplot!(:t, :T4, sampmarks=true, linealpha=0.4, nmarks=30)
@df thm_pd exptvwplot!(:t, :T3, nmarks=40, sampmarks=true, linealpha=0.4)
tendplot!(fitdat.t_end, label="", ls=:dash)
end
savefig(plotsdir("M1_exp.svg"))
