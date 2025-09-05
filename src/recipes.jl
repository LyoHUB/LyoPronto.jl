export qrf_integrate

# This is a plot recipe used to add markers only every so often along the series.
# e.g. if you have 100 data points but only want 10 markers, this will help.
# This is used inside many of the below plot recipes.
# Borrowed with some modification from https://github.com/JuliaPlots/Plots.jl/issues/2523#issuecomment-607090470
@recipe function f(::Type{Val{:samplemarkers}}, x, y, z; step = 10, offset=1)
    n = length(y)
    sx, sy = x[offset:step:n], y[offset:step:n]
    # add an empty series with the correct type for legend markers
    @series begin
        seriestype := :path
        markershape --> :auto
        x := [Inf]
        y := [Inf]
    end
    # add a series for the line
    @series begin
        primary := false # no legend entry
        markershape := :none # ensure no markers
        seriestype := :path
        seriescolor := get(plotattributes, :seriescolor, :auto)
        x := x
        y := y
    end
    # return  a series for the sampled markers
    primary := false
    seriestype := :scatter
    markershape --> :auto
    x := sx
    y := sy
end

@doc raw"""
    exptfplot(time, T1, [T2, ...]; nmarks=10, labsuffix=", exp.")
    exptfplot!(time, T1, [T2, ...]; nmarks=10, labsuffix=", exp.")

Plot recipe for one or more experimentally measured product temperatures, all at same times.
This recipe adds one series for each passed temperature series, with labels defaulting to `"T_{fi}"*labsuffix`.
"""
exptfplot
@doc (@doc exptfplot) exptfplot!

@userplot ExpTfPlot
@recipe function f(tpe::ExpTfPlot; nmarks=10, sampmarks=false, labsuffix = ", exp.")
    time, Ts... = tpe.args
    step = size(time, 1) ÷ nmarks
    n = size(Ts, 1)
    # pal = palette(:Blues_5,)[end:-1:begin] # Requires Plots for palette, so hard-code this default
    pal = [RGB{Float64}(0.031,0.318,0.612), 
           RGB{Float64}(0.192,0.51,0.741), 
           RGB{Float64}(0.42,0.682,0.839), 
           RGB{Float64}(0.741,0.843,0.906)]
    markers = (:circle, :square, :diamond, :hexagon)
    
    if length(Ts) == 1
        labels = ["\$T_\\mathrm{f}\$"*labsuffix]
    else
        labels = ["\$T_\\mathrm{f$i}\$"*labsuffix for i in 1:n]
    end
    for (i, T) in enumerate(Ts)
        @series begin
            markershape --> markers[i]
            label --> labels[i]
            seriescolor --> pal[i]
            if sampmarks
                seriestype := :samplemarkers
                step := step
                offset := step÷n *(i-1) + 1
                return time, T
            else
                seriestype --> :scatter
                offset = step÷n *(i-1) + 1
                minlen = min(length(time), length(T))
                return time[offset:step:minlen], T[offset:step:minlen]
            end
        end
    end
end

@doc raw"""
    exptvwplot(time, T1, [T2, ...]; trim=10, nmarks=10, labsuffix=", exp.")
    exptvwplot!(time, T1, [T2, ...]; trim=10, nmarks=10, labsuffix=", exp.")

Plot recipe for a set of experimentally measured vial wall temperatures.
This recipe adds one series for each passed temperature series, with labels defaulting to `"T_{vwi}"*labsuffix`.
`trim` is an integer, indicating how many points to skip at a time, so that 
the dotted line looks dotted even with noisy data.
"""
exptvwplot
@doc (@doc exptvwplot) exptvwplot!

@userplot ExpTvwPlot
@recipe function f(tpev::ExpTvwPlot; nmarks=30, sampmarks=false, trim=1, labsuffix=", exp.")
    time, Ts... = tpev.args
    n = size(Ts, 1)
    # color = palette(:Blues_5)[end]
    pal = [RGB{Float64}(0.031,0.318,0.612), 
           RGB{Float64}(0.192,0.51,0.741), 
           RGB{Float64}(0.42,0.682,0.839), 
           RGB{Float64}(0.741,0.843,0.906)]
    # markers = (:pentagon, :ltriangle, :rtriangle, :heptagon)
    markers = (:circle, :square, :diamond, :hexagon)
    
    if length(Ts) == 1
        labels = ["\$T_\\mathrm{vw}\$"*labsuffix]
    else
        labels = ["\$T_\\mathrm{vw$i}\$"*labsuffix for i in 1:n]
    end
    for (i, T) in enumerate(Ts)
        @series begin
            markershape --> markers[i]
            label --> labels[i]
            seriescolor --> pal[i]
            markercolor --> :white
            markerstrokecolor --> pal[i]
            markerstrokewidth --> 1.5
            if sampmarks
                linestyle --> :dash
                time_trim = time[begin:trim:end]
                T_trim = T[begin:trim:end]
                seriestype := :samplemarkers
                step = size(time_trim, 1) ÷ nmarks
                step := step
                offset := step÷n *(i-1) + 1
                return time_trim, T_trim
            else
                minlen = min(length(time), length(T))
                seriestype --> :scatter
                step = size(time, 1) ÷ nmarks
                offset = step÷n *(i-1) + 1
                return time[offset:step:minlen], T[offset:step:minlen]
            end
        end
    end
    
end

@doc raw"""
    modconvtplot(sols; sampmarks=false, labsuffix = ", model")
    modconvtplot!(sols; sampmarks=false, labsuffix = ", model")

Plot recipe for one or multiple solutions to the Pikal model, e.g. the output of [`gen_sol_pd`](@ref LyoPronto.gen_sol_pd).
This adds a series to the plot for each passed solution, with labels defaulting to `"T_{fi}"*labsuffix`.

If `sampmarks` is true, the solution will be interpolated to evenly spaced time points and 
markers will be added to some of those points.
"""
modconvtplot
@doc (@doc modconvtplot) modconvtplot!


@userplot ModConvTPlot
@recipe function f(tpmc::ModConvTPlot; sampmarks=false, labsuffix = ", model")
    sols = tpmc.args
    # pal = palette(:Oranges_4).colors[end:-1:begin+1] # Requires Plots as dependency...
    pal = [
    RGB{Float64}(0.851,0.278,0.004)
    RGB{Float64}(0.992,0.553,0.235)
    RGB{Float64}(0.992,0.745,0.522)
    ]
    markers = (:utriangle, :dtriangle, :cross)
    if length(sols) == 1
        labels = ["\$T_\\mathrm{f}\$"*labsuffix]
    else
        labels = ["\$T_\\mathrm{f$i}\$"*labsuffix for i in 1:length(sols)]
    end

    if sampmarks 
        for (i, sol) in enumerate(sols)
            t_nd = range(sol.t[begin], sol.t[end], length=101)
            time = t_nd*u"hr"
            T = sol.(t_nd, idxs=2)*u"K"
            @series begin
                seriestype := :samplemarkers
                step := 20
                markershape --> :auto
                seriescolor --> pal[i]
                label --> labels[i]
                return time, T
            end
        end
    else
        for (i, sol) in enumerate(sols)
            t_nd = sol.t
            time = t_nd*u"hr"
            T = sol[2,:]*u"K"
            @series begin
                seriescolor --> pal[i]
                label --> labels[i]
                return time, T
            end
        end
    end
end


"""
    modrftplot(sol, labsuffix=", model", sampmarks=false, trimend=0)
    modrftplot!(sol, labsuffix=", model", sampmarks=false, trimend=0)

Plot recipe for one solution to the lumped capacitance model.

This adds two series to the plot, with labels defaulting to `["T_f" "T_{vw}"] .* labsuffix`.
The optional argument `trimend` controls how many time points to trim from the end
(which is helpful if temperature shoots up as mf -> 0).

If `sampmarks` is true, the solution will be interpolated to evenly spaced time points and 
markers will be added to some of those points.

Since this is a recipe, any Plots.jl keyword arguments can be passed to modify the plot.
"""
modrftplot
@doc (@doc modrftplot) modrftplot!

@userplot ModRFTPlot
@recipe function f(tpmr::ModRFTPlot; labsuffix = ", model", sampmarks=false, trimend=0)
    sol = tpmr.args[1]
    if sampmarks
        t_nd = range(sol.t[begin], sol.t[end-trimend], length=31)
        step = 6
    else
        t_nd = sol.t[1:end-trimend]
    end
    time = t_nd*u"hr"
    # color = palette(:Oranges_3)[end]
    color = RGB{Float64}(0.902,0.333,0.051)
    # Frozen temperature: tends to have a crazy time point at end
    if sampmarks
        @series begin
            seriestype := :samplemarkers
            step := step
            T = sol.(t_nd, idxs=2)*u"K"
            markershape --> :dtriangle
            seriescolor --> color
            label --> "\$T_\\mathrm{f}\$"*labsuffix # Default label
            time, T
        end
        @series begin
            seriestype := :samplemarkers
            step := step
            T = sol.(t_nd, idxs=3)*u"K"
            markershape --> :dtriangle
            seriescolor --> color
            markercolor --> :white
            markerstrokecolor --> color
            label --> "\$T_\\mathrm{vw}\$"*labsuffix # Default label
            linestyle := :dash
            time, T
        end
    else
        @series begin
            T = sol.(t_nd, idxs=2)*u"K"
            seriescolor --> color
            label --> "\$T_\\mathrm{f}\$"*labsuffix # Default label
            time, T
        end
        @series begin
            T = sol.(t_nd, idxs=3)*u"K"
            seriescolor --> color
            label --> "\$T_\\mathrm{vw}\$"*labsuffix # Default label
            linestyle := :dash
            time, T
        end
    end
end

@doc raw"""
    tendplot(t_end)
    tendplot!(t_end)
    tendplot(t_end1, t_end2)
    tendplot!(t_end1, t_end2)

Plot recipe that adds a labeled vertical line to the plot at time `t_end`. 
A default label and styling are applied, but these can be modified by keyword arguments as usual
If two time points are passed, a light shading is applied between instead of a vertical line.
"""
tendplot
@doc (@doc tendplot) tendplot!

@userplot tendPlot
@recipe function f(tp::tendPlot)
    if ismissing(tp.args[1])
        return nothing
    elseif length(tp.args) == 1 && ~(tp.args[1] isa Tuple)
        t_end = tp.args[1]
        @series begin
            seriestype := :vline
            label --> "\$t_\\mathrm{end}\$"
            seriescolor --> :gray
            linestyle --> :dot
            linewidth --> 3
            return [ustrip(u"hr", t_end)]
        end
    elseif length(tp.args) == 2 || tp.args[1] isa Tuple
        if tp.args[1] isa Tuple
            t_end = tp.args[1]
        else
            t_end = tp.args[1:2]
        end
        @series begin
            seriestype := :vspan
            label --> "\$t_\\mathrm{end}\$, range"
            seriescolor --> :gray
            seriesalpha --> 0.3
            linewidth --> 1
            return [ustrip.(u"hr", t_end)...]
        end
    else
        error("tendPlot requires 1 or 2 arguments")
    end
end

@recipe function f(rv::RampedVariable{true, T1,T2,T3,T4}; tmax = rv.timestops[end]*100) where {T1,T2,T3,T4}
    t = vcat(rv.timestops, tmax)
    v = rv.(t)
    @series begin
        seriestype := :path
        return t, v
    end
end
@recipe function f(rv::RampedVariable{false, T1,T2,T3,T4}) where {T1,T2,T3,T4}
    @series begin
        seriestype := :hline
        return [rv(0)]
    end
end

@recipe function f(pdf::PrimaryDryFit)
    @series begin
        return ExpTfPlot((pdf.t, pdf.Tfs...))
    end
    if !ismissing(pdf.Tvws)
        if pdf.Tvws isa Number
            @series begin
                seriestype := :scatter
                label --> "\$T_\\mathrm{vw}\$"
                return [pdf.t[maximum(pdf.Tf_iend)]], [pdf.Tvws]
            end
        else 
            @series begin
                return ExpTvwPlot((pdf.t,pdf.Tvws...))
            end
        end
    end
    if !ismissing(pdf.t_end)
        @series begin
            return tendPlot(pdf.t_end)
        end
    end
end

@userplot AreaStackPlot
@recipe function f(asp::AreaStackPlot)
    x, ys... = asp.args
    yprev = zeros(eltype(ys[1]), size(ys[1], 1))
    ymat = hcat(yprev, ys...)
    ys_cum = cumsum(ymat, dims=2)
    pal = [:yellow, :orange, :red]
    for i in axes(ys_cum, 2)[begin:end-1]
        @series begin
            fillrange := ys_cum[:,i]
            fillcolor --> pal[i]
            linealpha --> 0
            marker --> :none
            return x,  ys_cum[:,i+1]
        end
        nothing
    end
end

@userplot BarStackPlot
@recipe function f(bsp::BarStackPlot)
    x, ys... = bsp.args
    yprev = zeros(eltype(ys[1]), size(ys[1], 1))
    ymat = hcat(yprev, ys...)
    ys_cum = cumsum(ymat, dims=2)
    ys_cum = ys_cum ./ maximum(ys_cum, dims=2)
    pal = [:yellow, :orange, :red]
    for i in axes(ys_cum, 2)[begin:end-1]
        @series begin
            seriestype := :bar
            fillrange := ys_cum[:,i]
            fillcolor --> pal[i]
            linewidth --> 1
            linecolor --> :black
            marker --> :none
            return x,  ys_cum[:,i+1]
        end
        nothing
    end
end

@doc raw"""
    qrf_integrate(sol, RF_params)

Compute the integral over time of each heat transfer mode in the lumped capacitance model.

Returns as a Dict{String, Quantity{...}}, with string keys `Qsub, Qshf, Qvwf, QRFf, QRFvw`.
"""
function qrf_integrate(sol, RF_params)

    # # This attempt currently falls down because it tries to materialize a zero vector and can't put zeros into a Unitful array
    # integrand = (u,t,integ) -> lumped_cap_rf(u, RF_params, t, energy_output=true)[2]
    # ex_res = integrand(sol.u0, 0, nothing)
    # integ_values = IntegrandValuesSum(ex_res)
    # integ_callback = IntegratingSumCallback(integrand, integ_values, ex_res.*0)
    # cbs = CallbackSet(integ_callback, end_drying_callback)
    # prob = ODEProblem(lumped_cap_rf, sol.u0, (0, 1e10), RF_params; callback=cbs)
    # sol = solve(prob, Rodas3(autodiff=false))

    # So we do a manual Riemann integration on the solution output
    Qcontrib = map(sol.t) do ti
        lumped_cap_rf!([0.0, 0.0, 0.0], sol(ti), RF_params, ti, Val(true))
    end
    Qcontrib = hcat(Qcontrib...)
    Qsub = Qcontrib[1,:]
    Qshf = Qcontrib[2,:]
    Qvwf = Qcontrib[3,:]
    QRFf = Qcontrib[4,:]
    QRFvw = Qcontrib[5,:]

    t = sol.t*u"hr"
    weights = fill(t[1], length(sol.t))
    dt = (t[begin+1:end] .- t[begin:end-1])
    weights[begin:end-1] += dt./2
    weights[begin+1:end] += dt./2

    qinteg = map([Qsub, Qshf, Qvwf, QRFf, QRFvw]) do q
        sum(q .* weights)
    end
    return Dict("Qsub"=>qinteg[1], 
                "Qshf"=>qinteg[2],
                "Qvwf"=>qinteg[3],
                "QRFf"=>qinteg[4],
                "QRFvw"=>qinteg[5])
end

"""
    pressurenames(name)

Attempt to catch common names or shorthands for Pirani and capacitance manometer, then give back a nice label.
"""
function pressurenames(name::String)
    if ~isnothing(match(r"[pP]ir", name))
        return "\$p_\\mathrm{ch}, \\mathrm{Pirani}\$"
    elseif !isnothing(match(r"[cC][mM]|[cC]ap", name)) 
        return "\$p_\\mathrm{ch}, \\mathrm{CM}\$"
    else
        return "\$p_\\mathrm{ch}\$"
    end
end
pressurenames(name) = pressurenames(string(name))

"""
    exppplot(t, p1, [p2, ...,], (name1, [name2, ...]))
    exppplot!(t, p1, [p2, ...,], (name1, [name2, ...]))

Plot recipe for plotting pressures of a lyophilization cycle.
"""
exppplot
@doc (@doc exppplot) exppplot!

@userplot ExpPPlot
@recipe function f(epp::ExpPPlot)
    t, ps..., names = epp.args
    defaultlabels = [pressurenames(name) for name in names]
    for (i, p) in enumerate(ps)
        @series begin
            seriestype --> :samplemarkers
            step --> length(p) ÷ 10
            ylabel --> "Pressure"
            label --> defaultlabels[i]
            return t, p
        end
    end
end

"""
    cycledataplot(tablelike, (Tname1, ...), Tsh_name, (p_name1, ...))
    cycledataplot!(tablelike, (Tname1, ...), Tsh_name, (p_name1, ...))

Plot recipe for plotting both temperatures and pressures of a lyophilization cycle.

This requires at least two subplots, one for temperature and one for pressure.
Either call `twinx(plot())` then `cycledataplot!(...)`, or give a layout like 
`cycledataplot(..., layout=(2,1))`"
"""
cycledataplot
@doc (@doc cycledataplot) cycledataplot!

@userplot CycleDataPlot
@recipe function f(cdp::CycleDataPlot; pcolor=:black, tendkind=:onoff)

    dat, Ts, Tsh, ps = cdp.args
    unitformat --> :square
    if !isnothing(Tsh)
        @series begin
            lw --> 2
            subplot := 1
            seriescolor --> :black
            label --> "\$T_\\mathrm{sh}\$"
            dat.t, getproperty(dat, Tsh)
        end
    end
    @series begin
        subplot := 1
        ylabel --> "Temperature"
        xlabel --> "Time"
        ygrid --> true
        lw --> 2
        ExpTfPlot((dat.t, map(x->getproperty(dat, x), Ts)...))
    end
    if length(plotattributes[:plot_object].subplots) < 2
        @warn "CycleDataPlot requires at least two subplots, one for temperature and one for pressure.
        Either call `twinx(plot())` then `cycledataplot!(...)`, 
        or give a layout like `cycledataplot(..., layout=(2,1))`"
    end
    subplot := 2
    foreground_color_axis --> pcolor
    bordercolor --> pcolor
    seriescolor --> pcolor
    @series begin
        ExpPPlot((dat.t, map(x->getproperty(dat, x), ps)..., ps))
    end
end