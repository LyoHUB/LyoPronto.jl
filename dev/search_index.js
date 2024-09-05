var documenterSearchIndex = {"docs":
[{"location":"alldocstrings/#All-Docstrings","page":"Reference","title":"All Docstrings","text":"","category":"section"},{"location":"alldocstrings/","page":"Reference","title":"Reference","text":"CurrentModule = LyoPronto","category":"page"},{"location":"alldocstrings/#Types","page":"Reference","title":"Types","text":"","category":"section"},{"location":"alldocstrings/","page":"Reference","title":"Reference","text":"Modules = [LyoPronto]\nPages = [\"structs.jl\"]","category":"page"},{"location":"alldocstrings/#LyoPronto.RampedVariable","page":"Reference","title":"LyoPronto.RampedVariable","text":"A convenience type for computing temperatures, pressures, etc. with multiple setpoints in sequence,  and linear interpolation according to a fixed ramp rate between set points\n\nThree main constructors are available: For a non-varying value, call with one argument: RampedVariable(constant_setpt) For one ramp from initial value to set point with indefinite hold, call with two arguments: RampedVariable(setpts, ramprate) And for multiple setpoints, call with three arguments: RampedVariable(setpts, ramprates, holds)\n\nWith three arguments, setpts, ramprates, and holds should all be vectors, with lengths N+1, N, N-1 respectively.\n\nThe resulting RampedVariable rv = RampedVariable(...) can be called as rv(x) at any value of x,  and will return the value at that time point along the ramp process.\n\n\n\n\n\n","category":"type"},{"location":"alldocstrings/#LyoPronto.RpFormFit","page":"Reference","title":"LyoPronto.RpFormFit","text":"A convenience type for dealing with the common functional form given to Rp and Kv.\n\nAn object Rp = RpFormFit(A, B, C) can be called as Rp(x), which simply computes A + B*x/(1 + C*x). Likewise, Kv = RpFormFit(Kc, Kp, Kd) can be called as Kv(p) to get Kc + Kp*p/(1 + Kd*p).\n\n\n\n\n\n","category":"type"},{"location":"alldocstrings/#Parameter-Fitting","page":"Reference","title":"Parameter Fitting","text":"","category":"section"},{"location":"alldocstrings/","page":"Reference","title":"Reference","text":"Modules = [LyoPronto]\nPages = [\"paramfits.jl\"]","category":"page"},{"location":"alldocstrings/#LyoPronto.gen_sol_conv_dim-NTuple{4, Any}","page":"Reference","title":"LyoPronto.gen_sol_conv_dim","text":"gen_sol_conv_dim(KRp_prm, otherparams, u0, tspan; kwargs...)\n\nSolve the Pikal model for primary drying with given Kv & Rp, returning the solution object and the set of parameters passed to solve.\n\notherparams contains: (hf0, c_solid, ρ_solution, Av, Ap, pch, Tsh)\nu0 is given as floats (not Unitful quantities), with dimensions [cm, K].\nKRp_prm has the form [Kv, R0, A1, A2]; this function assigns the units [W/m^2/K, cm^2Torrhr/g, cmTorrhr/g, 1/cm] to those numbers\nkwargs is passed directly (as is) to the ODE solve call.\n\nThis is used as a helper for obj_tT_conv to assemble a function which takes [Kv, R0, A1, A2] and returns sum squared error against experiment. To that end, use this to make an anonymous function of the form gen_sol = x->gen_sol_conv_dim(x, case_params, case_u0, case_tspan). That function is what you will pass to obj_tT_conv.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.gen_sol_rf_dim-NTuple{4, Any}","page":"Reference","title":"LyoPronto.gen_sol_rf_dim","text":"gen_sol_rf_dim(fitprm, params_bunch, u0, tspan; kwargs...)\n\nSolve the lumped-capacitance model for microwave-assisted primary drying with given fit parameters, returning the solution object and the set of parameters passed to solve.\n\nfitprm has the form [α, Kvwf, Bf, Bvw]; this function assigns the units [cm^1.5, cal/s/K/cm^2, Ω/m^2, Ω/m^2] to those numbers.\nparams_bunch contains a full listing of parameters used for the model, according to lumped_cap_rf, including dummy values for the fit parameters.\nu0 is given as floats (not Unitful quantities), with dimensions [g, K, K].\nkwargs is passed directly (as is) to the ODE solve call.\n\nThis is used as a helper for obj_tT_rf, obj_tTT_rf, obj_ttTT_rf to assemble a function which takes [α, Kvwf, Bf, Bvw] and returns sum squared error against experiment. To that end, use this to make an anonymous function of the form gen_sol = x->gen_sol_rf_dim(x, case_params, case_u0, case_tspan). That function is what you will pass to obj_tT_rf (or similar).\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.obj_tTT_rf-Tuple{Any, Any, Any}","page":"Reference","title":"LyoPronto.obj_tTT_rf","text":"obj_tTT_rf(fitprm, gen_sol, tTTdat; t_end=0.0u\"hr\", tweight=1, verbose=true)\n\nEvaluate an objective function which compares model solution with fitprm to experimental data in tTTdat.\n\ngen_sol is a function taking [α, Kvwf, Bf, Bvw] and returning a solution to lumped-capacitance microwave-assisted model; see gen_sol_rf_dim.\ntTTdat is experimental temperature series, of the form (time, Tf, Tvw), so with frozen and vial wall temperatures taken at the same time points.    See also obj_tT_rf and obj_ttTT_rf.\ntweight gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.\nt_end has a default value of 0.0u\"hr\", which (if left at default) is replaced with the last given time point.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.obj_tT_conv-Tuple{Any, Any, Any}","page":"Reference","title":"LyoPronto.obj_tT_conv","text":"obj_tT_conv(KRp_prm, gen_sol, tTdat; t_end=0.0u\"hr\", tweight=1, verbose = true)\n\nEvaluate an objective function which compares model solution with KRp_prm to experimental data in tTdat.\n\nArguments:\n\ngen_sol is a function taking [Kv, R0, A1, A2] and returning a solution to Pikal model; see gen_sol_conv_dim.\ntTdat is experimental temperature series, of the form (time, Tf).\ntweight gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.\nt_end has a default value of 0.0u\"hr\", which (if left at default) is replaced with the last time point.\nverbose defaults to true, in which case each call to this function prints info on the passed parameters, etc.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.obj_tT_rf-Tuple{Any, Any, Any}","page":"Reference","title":"LyoPronto.obj_tT_rf","text":"obj_tT_rf(fitprm, gen_sol, tTdat; t_end=0.0u\"hr\", tweight=1, Tvw_end = 0.0u\"K\", verbose=true)\n\nEvaluate an objective function which compares model solution with fitprm to experimental data in tTdat.\n\ngen_sol is a function taking [α, Kvwf, Bf, Bvw] and returning a solution to lumped-capacitance microwave-assisted model; see gen_sol_rf_dim.\ntTTdat is experimental temperature series, of the form (time, Tf), so with frozen temperatures only.    See also obj_tTT_rf and obj_ttTT_rf.\nTvw_end defaults to 0.0u\"K\", in which case vial wall temperatures are excluded from the objective.    If another value is passed, then the final model vial wall temperature is compared to that value and included in the objective.\ntweight gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.\nt_end has a default value of 0.0u\"hr\", which (if left at default) is replaced with the last given time point.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.obj_ttTT_rf-Tuple{Any, Any, Any}","page":"Reference","title":"LyoPronto.obj_ttTT_rf","text":"obj_ttTT_rf(fitprm, gen_sol, tTTdat; t_end=0.0u\"hr\", tweight=1, verbose=true)\n\nEvaluate an objective function which compares model solution with fitprm to experimental data in tTdat.\n\ngen_sol is a function taking [α, Kvwf, Bf, Bvw] and returning a solution to lumped-capacitance microwave-assisted model; see gen_sol_rf_dim.\ntTTdat is experimental temperature series, of the form (time_Tf, time_vw, Tf, Tvw), so with Tf and Tvw having separate time points.    This is useful if there is an early temperature rise in Tf, but Tvw continues to be reliable, so the model can fit to as much of Tvw as reasonable.   See also obj_tTT_rf and obj_tT_rf.\ntweight gives the weighting (in K^2/hr^2) of the end of drying in the objective, as compared to the temperature error.\nt_end has a default value of 0.0u\"hr\", which (if left at default) is replaced with the last given time point.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#Plot-Recipes","page":"Reference","title":"Plot Recipes","text":"","category":"section"},{"location":"alldocstrings/","page":"Reference","title":"Reference","text":"Modules = [LyoPronto]\nPages = [\"recipes.jl\"]","category":"page"},{"location":"alldocstrings/#LyoPronto.qrf_integrate-Tuple{Any, Any}","page":"Reference","title":"LyoPronto.qrf_integrate","text":"qrf_integrate(sol, RF_params)\n\nCompute the integral over time of each heat transfer mode in the lumped capacitance model.\n\nReturns as a Dict{String, Quantity{...}}, with string keys Qsub, Qshf, Qvwf, QRFf, QRFvw.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.tendplot","page":"Reference","title":"LyoPronto.tendplot","text":"tendplot(t_end)\ntendplot!(t_end)\n\nPlot recipe that adds a labeled vertical line to the plot at time t_end.  A default label and styling are applied, but these can be modified by keyword arguments as usual\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#LyoPronto.tendplot!","page":"Reference","title":"LyoPronto.tendplot!","text":"tendplot(t_end)\ntendplot!(t_end)\n\nPlot recipe that adds a labeled vertical line to the plot at time t_end.  A default label and styling are applied, but these can be modified by keyword arguments as usual\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#LyoPronto.tplotexperimental","page":"Reference","title":"LyoPronto.tplotexperimental","text":"tplotexperimental(time, T1, [T2, ...])\ntplotexperimental!(time, T1, [T2, ...])\n\nPlot recipe for one or more experimentally measured product temperatures, all at same times. This recipe adds one series for each passed temperature series, so pass labels as appropriate.\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#LyoPronto.tplotexperimental!","page":"Reference","title":"LyoPronto.tplotexperimental!","text":"tplotexperimental(time, T1, [T2, ...])\ntplotexperimental!(time, T1, [T2, ...])\n\nPlot recipe for one or more experimentally measured product temperatures, all at same times. This recipe adds one series for each passed temperature series, so pass labels as appropriate.\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#LyoPronto.tplotexpvw","page":"Reference","title":"LyoPronto.tplotexpvw","text":"tplotexpvw(time, temperature)\ntplotexpvw!(time, temperature)\n\nPlot recipe for a set of experimentally measured vial wall temperatures. This recipe adds only one series to the plot.\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#LyoPronto.tplotexpvw!","page":"Reference","title":"LyoPronto.tplotexpvw!","text":"tplotexpvw(time, temperature)\ntplotexpvw!(time, temperature)\n\nPlot recipe for a set of experimentally measured vial wall temperatures. This recipe adds only one series to the plot.\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#LyoPronto.tplotmodelconv","page":"Reference","title":"LyoPronto.tplotmodelconv","text":"tplotmodelconv(sols)\ntplotmodelconv!(sols)\n\nPlot recipe for one or multiple solutions to the Pikal model, e.g. the output of gen_sol_conv_dim. This adds one series to the plot for each passed solution, so pass as many labels (e.g. [\"Tf1\" \"Tf2\"]) to this plot call as solutions to add labels to the legend.\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#LyoPronto.tplotmodelconv!","page":"Reference","title":"LyoPronto.tplotmodelconv!","text":"tplotmodelconv(sols)\ntplotmodelconv!(sols)\n\nPlot recipe for one or multiple solutions to the Pikal model, e.g. the output of gen_sol_conv_dim. This adds one series to the plot for each passed solution, so pass as many labels (e.g. [\"Tf1\" \"Tf2\"]) to this plot call as solutions to add labels to the legend.\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#LyoPronto.tplotmodelrf","page":"Reference","title":"LyoPronto.tplotmodelrf","text":"tplotmodelrf(sol)\ntplotmodelrf!(sol)\n\nPlot recipe for one solution to the lumped capacitance model, e.g. the output of gen_sol_rf_dim. This adds two series to the plot, so pass two labels (e.g. [\"Tf\" \"Tvw\"]) to this plot call to add labels to the legend.\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#LyoPronto.tplotmodelrf!","page":"Reference","title":"LyoPronto.tplotmodelrf!","text":"tplotmodelrf(sol)\ntplotmodelrf!(sol)\n\nPlot recipe for one solution to the lumped capacitance model, e.g. the output of gen_sol_rf_dim. This adds two series to the plot, so pass two labels (e.g. [\"Tf\" \"Tvw\"]) to this plot call to add labels to the legend.\n\n\n\n\n\n","category":"function"},{"location":"alldocstrings/#Vial-Dimensions","page":"Reference","title":"Vial Dimensions","text":"","category":"section"},{"location":"alldocstrings/","page":"Reference","title":"Reference","text":"Modules = [LyoPronto]\nPages = [\"get_vial_dims.jl\"]","category":"page"},{"location":"alldocstrings/#LyoPronto.get_vial_mass-Tuple{String}","page":"Reference","title":"LyoPronto.get_vial_mass","text":"get_vial_mass(vialsize::String)\n\nReturn vial mass for given ISO vial size.\n\nUses a table provided by Schott, stored internally in a CSV.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.get_vial_radii-Tuple{String}","page":"Reference","title":"LyoPronto.get_vial_radii","text":"get_vial_radii(vialsize::String)\n\nReturn inner and outer radius for passed ISO vial size.\n\nUses a table provided by Schott, stored internally in a CSV.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.get_vial_shape-Tuple{String}","page":"Reference","title":"LyoPronto.get_vial_shape","text":"get_vial_shape(vialsize::String)\n\nReturn a Dict{Symbol, Any} with a slew of vial dimensions, useful for drawing the shape of the vial with make_outlines.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.get_vial_thickness-Tuple{String}","page":"Reference","title":"LyoPronto.get_vial_thickness","text":"get_vial_thickness(vialsize::String)\n\nReturn vial wall thickness for given ISO vial size.\n\nUses a table provided by Schott, stored internally in a CSV.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.make_outlines-Tuple{Any, Any}","page":"Reference","title":"LyoPronto.make_outlines","text":"make_outlines(dims, Vfill)\n\nReturn a sequence of points (ready to be made into Plots.Shapes for the vial and fill volume, with Unitful dimensions, for given vial dimensions.\n\nThis is a convenience function for making figures illustrating fill depth.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#Model-Equations","page":"Reference","title":"Model Equations","text":"","category":"section"},{"location":"alldocstrings/","page":"Reference","title":"Reference","text":"Modules = [LyoPronto]\nPages = [\"model.jl\"]","category":"page"},{"location":"alldocstrings/#LyoPronto.end_drying_callback","page":"Reference","title":"LyoPronto.end_drying_callback","text":"A callback for use in simulating either the Pikal or RF model.\n\nTerminates the time integration when end_cond evaluates to true.\n\n\n\n\n\n","category":"constant"},{"location":"alldocstrings/#LyoPronto.calc_psub-Tuple{F} where F<:Number","page":"Reference","title":"LyoPronto.calc_psub","text":"calc_psub(T::F) where F<:Number\ncalc_psub(T::Q) where Q<:Quantity\n\nCompute pressure (in Pascals) of sublimation at temperature T in Kelvin.\n\nThis is essentially an Arrhenius fit, where we compute: psub = pref * exp(-ΔHsub*Mw / R T)\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.end_cond-Tuple{Any, Any, Any}","page":"Reference","title":"LyoPronto.end_cond","text":"end_cond(u, t, integ)\n\nCompute the end condition for primary drying (that mf or hf approaches zero).\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.lumped_cap_rf-Tuple{Any, Any, Any}","page":"Reference","title":"LyoPronto.lumped_cap_rf","text":"lumped_cap_rf(u, params, tn)\n\nCompute the right-hand-side function for the ODEs making up the lumped-capacitance microwave-assisted model.\n\nSpecifically, this is [dmf/dt, dTf/dt, dTvw/dt] given u = [mf, Tf, Tvw]. u is taken without units but assumed to have the units of [g, K, K] (which is internally added). tn is assumed to be in hours (internally added), so dudt is returned with units [g/hr, K/hr, K/hr] to be consistent.\n\nThe full set of necessary parameters is given in the form of a tuple-of-tuples:\n\nparams = (   \n    (Rp, h_f0, cSolid, ρ_solution),\n    (K_shf_f, A_v, A_p),\n    (pch, Tsh, P_per_vial),\n    (m_f0, cp_f, m_v, cp_v, A_rad),\n    (f_RF, epp_f, epp_vw),\n    (alpha, K_vwf, B_f, B_vw),\n)\n\nThese should all be Unitful quantities with appropriate dimensions, with some exceptions which are callables returning quantities. See RpFormFit and `RampedVariable for convenience types that can help with these cases.\n\nRp(x) with x a length returns mass transfer resistance (as a Unitful quantity)\nK_shf_f(p) with p a pressure returns heat transfer coefficient (as a Unitful quantity).\nTsh(t), pch(t), P_per_vial(t) return shelf temperature, chamber pressure, and microwave power respectively at time t.\n\nFor implementation details, see lumped_cap_rf_model.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.lumped_cap_rf_model-Tuple{Any, Any, Any}","page":"Reference","title":"LyoPronto.lumped_cap_rf_model","text":"lumped_cap_rf_model(u, params, tn)\n\nThis does the work for lumped_cap_rf, but returns dudt,  [Q_sub, Q_shf, Q_vwf, Q_RF_f, Q_RF_vw, Q_shw] with Q_... as Unitful quantities in watts.  The extra results are helpful in investigating the significance of the various heat transfer modes.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.lyo_1d_dae!-NTuple{4, Any}","page":"Reference","title":"LyoPronto.lyo_1d_dae!","text":"lyo_1d_dae!(du, u, params, t)\n\nInternal implementation of the Pikal model. See lyo_1d_dae_f for the wrapped version, which is more fully documented.\n\n\n\n\n\n","category":"method"},{"location":"alldocstrings/#LyoPronto.lyo_1d_dae_f","page":"Reference","title":"LyoPronto.lyo_1d_dae_f","text":"lyo_1d_dae_f = ODEFunction(lyo_1d_dae!, mass_matrix=lyo_1d_mm)\n\nCompute the right hand side function for the Pikal model.\n\nThe DAE system which is the Pikal model (1 ODE, one nonlinear algebraic equation for pseudosteady conditions) is here treated as a constant-mass-matrix implicit ODE system. The implementation is in lyo_1d_dae!\n\nThe initial conditions u0 = [h_f, Tf] should be unitless, but are internally assigned to be in [cm, K]. The unitless time is taken to be in hours, so derivatives are given in unitless [cm/hr, K/hr].\n\nparams has the form:\n\nparams = (\n    (Rp, hf0, c_solid, ρ_solution),\n    (Kshf, Av, Ap),\n    (pch, Tsh) ,\n)\n\nwhere some quantities with Unitful units and some are callables returning quantities See RpFormFit and RampedVariable for convenience types that can help with the callables.\n\nRp(x) with x a length returns mass transfer resistance (as a Unitful quantity)\nK_shf_f(p) with p a pressure returns heat transfer coefficient (as a Unitful quantity).\nTsh(t), pch(t) return shelf temperature and chamber pressure respectively at time t.\n\n\n\n\n\n","category":"function"},{"location":"#Home","page":"Home","title":"Home","text":"","category":"section"},{"location":"#LyoPronto.jl","page":"Home","title":"LyoPronto.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Julia package providing common computations for pharmaceutical lyophilization.","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This relatively small package puts together some convenience functions useful for simulating primary drying in lyophilization.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"As a Julia package, this code can be easily installed with the Julia package manager. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"From the Julia REPL's Pkg mode (open a REPL and type ] so that the prompt turns blue), add this package as a Git repo:","category":"page"},{"location":"","page":"Home","title":"Home","text":"add https://github.com/LyoHUB/LyoPronto.jl.git","category":"page"},{"location":"","page":"Home","title":"Home","text":"dev can be substituted for add if you want to make changes to this package yourself, as explained in the Julia Pkg manual.","category":"page"},{"location":"#Dependencies-and-Reexports","page":"Home","title":"Dependencies and Reexports","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package leverages the strengths of the DifferentialEquations.jl ecosystem to solve equations quickly and efficiently, although it only directly depends on OrdinaryDiffEqRosenbrock and DiffEqCallbacks, which are both reexported.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Also provided are plot recipes for Plots.jl, although this package only depends on RecipesBase.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Heavy use is made of Unitful.jl, which is reexported.","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Written by Isaac S. Wheeler, a PhD student at Purdue University. This work was supported in part by funding for NIIMBL project PC4.1-307 .","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"None yet. My intentions are to use the MIT license once this has been published in a scientific journal.","category":"page"},{"location":"example/#Example:-Tuning-Mass-Transfer-Resistance","page":"Example","title":"Example: Tuning Mass Transfer Resistance","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"The source code for this example can be found in scripts/tuning_sucrose_Rp.jl. ","category":"page"},{"location":"example/#Setup","page":"Example","title":"Setup","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"My current recommended approach for making use of this package is to have installed LyoPronto as a dependency of your project, according to the install instructions","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"For management of your project, DrWatson.jl is an effective tool which I like. I will demonstrate one use of that tool here.","category":"page"},{"location":"example/#Packages-Used","page":"Example","title":"Packages Used","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"The following packages are used in this example:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"First, load packages:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"using DrWatson\nusing LyoPronto\nusing CSV\nusing Plots\nusing LaTeXStrings\nusing SavitzkyGolay\nusing Optim","category":"page"},{"location":"example/#Loading-Experimental-Data","page":"Example","title":"Loading Experimental Data","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"First, load a .csv of experimental data with the CSV.File method, and use DrWatson's datadir function to help keep things clean.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"This file has 7 rows of information at the top, so the data headers are at row 8; it is located in the project directory under [project]/data/exp_raw/\"2024-06-21-16_MFD_AH.csv.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"procdat = CSV.File(datadir(\"exp_raw\", \"2024-06-21-16_MFD_AH.csv\"), header=8)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"For the Millrock MicroFD lyophilizer that data was gathered on, the .csv file has a column Phase indicating what step of the process is going on, so to get our primary drying data out we can identify it by rows where the Phase == 4:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"is_PD = procdat[\"Phase\"] .== 4","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"To get the data series, we index to those rows for each variable of interest, and mark their units from Unitful.jl by multiplying by u\"[unit]\". The time points are once a minute, so here we construct a set of time data from scratch rather than mess with reading in the DateTime values from the lyophilizer (which may flip to zero at midnight and have other annoying behavior).","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"tproc = range(0, length=sum(is_PD), step=1/60)*u\"hr\"\n\np_pir = procdat[\"VacPirani\"][is_PD]*u\"mTorr\"\nTsh_d = procdat[\"ShelfSetPT\"][is_PD]*u\"°C\"\nT1 = procdat[\"TP1\"][is_PD]*u\"°C\"\nT2 = procdat[\"TP2\"][is_PD]*u\"°C\"\nT3 = procdat[\"TP4\"][is_PD]*u\"°C\"","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"To check that everything looks right, plot the temperatures, taking advantage of a recipe from this package, as well as the L\"[latex]\" macro from LaTeXStrings.  (This requires Plots.jl.)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"tplotexperimental(tproc, T1, T2, T3, labels=[L\"T_{p1}\" L\"T_{p2}\" L\"T_{p3}\"])\nplot!(tproc, Tsh_d, c=:black)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Based on an examination of the temperature data, we want to go only up to the \"temperature rise\" commonly observed in lyophilization near (but not at) the end of drying. This happens at 7.5 hours for T1 and T3 and at about 12.5 hours for T2.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"T1trm = T1[tproc .< 7.5u\"hr\"]\nT2trm = T2[tproc .< 12.5u\"hr\"]\nT3trm = T3[tproc .< 7.5u\"hr\"]","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Look at the Pirani pressure and ascertain the end of drying by using a Savitzky-Golay filter to identify a maximum in the second derivative. Note that because SavitzkyGolay doesn't play nice with Unitful, we strip out units and add them back in. Another plot recipe plots the end of drying as a vertical line on that pressure graph.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"p_pir_sm = savitzky_golay(ustrip.(u\"mTorr\", p_pir), 91, 3, deriv=0).y *u\"mTorr\" # 91 a window width; 3 the polynomial order\np_pir_der2 = savitzky_golay(ustrip.(u\"mTorr\", p_pir), 91, 3, deriv=2).y\nt_end = tproc[argmax(p_pir_der2[100:end])+99]\nplot(tproc, p_pir, label=\"data\")\nplot!(tproc, p_pir_sm, label=\"smoothed data\")\ntendplot!(t_end, label=L\"t_{end}\")","category":"page"},{"location":"example/#Set-up-model-parameters","page":"Example","title":"Set up model parameters","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"Below, we make liberal use of Unitful units. Also note that RampedVariable and RpFormFit are used to simplify some common things we use.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"# Vial geometry\nvialsize = \"6R\"\nAp, Av = π.*get_vial_radii(vialsize).^2\n\n# Formulation parameters; Rp here is a placeholder guess\nc_solid = 0.05u\"g/mL\" # g solute / mL solution\nρ_solution = 1u\"g/mL\" # g/mL total solution density\nR0 = 0.8u\"cm^2*Torr*hr/g\"\nA1 = 14u\"cm*Torr*hr/g\"\nA2 = 1u\"1/cm\"\nRp = RpFormFit(R0, A1, A2)\n\n# Fill\nVfill = 3u\"mL\" # ml\nhf0 = Vfill / Ap\n\n# Cycle parameters\npch = RampedVariable(70u\"mTorr\") # constant pressure\nT_shelf_0 = -40.0u\"°C\" # initial shelf temperature\nT_shelf_final = -10.0u\"°C\"  # final shelf temperature\nramp_rate = 0.5 *u\"K/minute\" # ramp rate\n# Ramp for shelf temperature: convert to Kelvin because Celsius doesn't do math very well\nTsh = RampedVariable(uconvert.(u\"K\", [T_shelf_0, T_shelf_final]), ramp_rate)\n\n# Guess for heat transfer\nKC = 6.556e-5u\"cal/s/K/cm^2\"\nKP = 2.41e-3u\"cal/s/K/cm^2/Torr\"\nKD = 2.62u\"1/Torr\"\nKshf = RpFormFit(KC, KP, KD)\n# Alternative guess, without pressure dependence:\n# Kshf = p-> 5.0u\"W/m^2/K\"\n\nparams_bunch = [\n    (Rp, hf0, c_solid, ρ_solution),\n    (Kshf, Av, Ap),\n    (pch, Tsh)\n]\n","category":"page"},{"location":"example/#Run-a-check-simulation","page":"Example","title":"Run a check simulation","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"# Time span: used to set initial time and to give an upper bound on time, in case parameters are bad\ntspan = (0.0, 100.0) # hours\n# Initial condition\nu0 = [ustrip(u\"cm\", hf0), ustrip(u\"K\", Tsh(0u\"minute\"))]\n\n# Set up as an ODE problem\nprob = ODEProblem(lyo_1d_dae_f, u0, tspan, tuple(params_bunch...))\n# Solve with the Rodas4P() algorithm, use a callback to terminate at end of drying\nsol = solve(prob, Rodas4P(), callback=LyoPronto.end_drying_callback)\n","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"And then plot these results with a recipe, to make sure that the chosen K_sh-f and R_p are sane starting points:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(tproc, Tsh_d, c=:black, label=L\"T_{sh}\")\ntplotmodelconv!(sol, label=L\"$T_p$, model\")","category":"page"},{"location":"example/#Minimize-least-square-difference-to-compare-to-experimental-data","page":"Example","title":"Minimize least square difference to compare to experimental data","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"First, set up a function that takes K_sh-f and R_p and returns a solution object (see gen_sol_conv_dim ).","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"otherparams = (hf0, c_solid, ρ_solution, Av, Ap, pch, Tsh)\ngen_sol_conv = KRp -> gen_sol_conv_dim(KRp, otherparams, u0, tspan)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Next, set up an objective function we will minimize (see obj_tT_conv for details). We will compare to T1.","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"obj_KRp1 = KRp -> obj_tT_conv(KRp, gen_sol_conv, (tproc[tproc.<7.5u\"hr\"], T1trm), t_end=t_end) ","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Minimize this objective function, using Optim with the NelderMead() algorithm:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"opt_KRp1 = optimize(obj_KRp1, [12, 0.1, 5, 0.1], NelderMead())","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Get out the found values of our tuning parameters, and generate the corresponding solution profile:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"KRp_1 = Optim.minimizer(opt_Krp1)\nprof1 = gen_sol_conv(Krp_1)[1]","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"Plot and compare experiment to model:","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"tplotexperimental(tproc, T1, T2, T3)\ntplotmodelconv!(prof1)\ntendplot!(t_end)","category":"page"}]
}
