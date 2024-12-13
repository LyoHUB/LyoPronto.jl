export qrf_integrate

# This is a plot recipe used to add markers only every so often along the series.
# e.g. if you have 100 data points but only want 10 markers, this will help.
# This is used inside many of the below plot recipes.
@recipe function f(::Type{Val{:samplemarkers}}, x, y, z; step = 10, offset=1)
    n = length(y)
    sx, sy = x[offset:step:n], y[offset:step:n]
    linewidth --> 2.5
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
    exptfplot(time, T1, [T2, ...]; labsuffix=", exp.")
    exptfplot!(time, T1, [T2, ...]; labsuffix=", exp.")

Plot recipe for one or more experimentally measured product temperatures, all at same times.
This recipe adds one series for each passed temperature series, with labels defaulting to `"T_{fi}"*labsuffix`.
"""
exptfplot
@doc (@doc exptfplot) exptfplot!

@userplot ExpTfPlot
@recipe function f(tpe::ExpTfPlot; labsuffix = ", exp.")
    time, Ts... = tpe.args
    step = size(time, 1) ÷ 10
    n = size(Ts, 1)
    # pal = palette(:Blues_5,)[end:-1:begin] # Requires Plots for palette, so hard-code this default
    pal = [RGB{Float64}(0.031,0.318,0.612), 
           RGB{Float64}(0.192,0.51,0.741), 
           RGB{Float64}(0.42,0.682,0.839), 
           RGB{Float64}(0.741,0.843,0.906)]
    
    if length(Ts) == 1
        labels = ["\$T_{f}\$"*labsuffix]
    else
        labels = ["\$T_{f$i}\$"*labsuffix for i in 1:n]
    end
    for (i, T) in enumerate(Ts)
        minlen = min(length(time), length(T))
        @series begin
            seriestype := :samplemarkers
            step := step
            offset := step÷n *(i-1) + 1
            # markershape --> :auto
            markersize --> 7
            label --> labels[i]
            seriescolor --> pal[i]
            if T == :dummy
                return [Inf], [Inf]
            else
                return time[1:minlen], T[1:minlen]
            end
        end
    end
end

@doc raw"""
    exptvwplot(time, T1, [T2, ...]; trim)
    exptvwplot!(time, T1, [T2, ...]; trim)

Plot recipe for a set of experimentally measured vial wall temperatures.
This recipe adds one series for each passed temperature series, with labels defaulting to `"T_{vwi}"*labsuffix`.
`trim` is an integer, indicating how many points to skip at a time, so that 
the dotted line looks dotted even with noisy data.
"""
exptvwplot
@doc (@doc exptvwplot) exptvwplot!

@userplot ExpTvwPlot
@recipe function f(tpev::ExpTvwPlot; trim=10, labsuffix=", exp.")
    time, Ts... = tpev.args
    time_trim = time[begin:trim:end]
    step = size(time_trim, 1) ÷ (10)
    n = size(Ts, 1)
    # color = palette(:Blues_5)[end]
    pal = [RGB{Float64}(0.031,0.318,0.612), 
           RGB{Float64}(0.192,0.51,0.741), 
           RGB{Float64}(0.42,0.682,0.839), 
           RGB{Float64}(0.741,0.843,0.906)]
    
    if length(Ts) == 1
        labels = ["\$T_{vw}\$"*labsuffix]
    else
        labels = ["\$T_{vw$i}\$"*labsuffix for i in 1:n]
    end
    for (i, T) in enumerate(Ts)
        T_trim = T[begin:trim:end]
        minlen = min(length(time_trim), length(T_trim))
        @series begin
            seriestype := :samplemarkers
            step := step
            offset := step÷n *(i-1) + 1
            linestyle --> :dash
            markersize --> 7
            label --> labels[i]
            seriescolor --> pal[i]
            return time_trim[1:minlen], T_trim[1:minlen]
        end
    end
    
end

@doc raw"""
    modconvtplot(sols; labsuffix = ", model")
    modconvtplot!(sols; labsuffix = ", model")

Plot recipe for one or multiple solutions to the Pikal model, e.g. the output of [`gen_sol_conv_dim`](@ref LyoPronto.gen_sol_conv_dim).
This adds a series to the plot for each passed solution, with labels defaulting to `"T_{fi}"*labsuffix`.
"""
modconvtplot
@doc (@doc modconvtplot) modconvtplot!


@userplot ModConvTPlot
@recipe function f(tpmc::ModConvTPlot; labsuffix = ", model")
    sols = tpmc.args
    # pal = palette(:Oranges_4).colors[end:-1:begin+1] # Requires Plots as dependency...
    pal = [
    RGB{Float64}(0.851,0.278,0.004)
    RGB{Float64}(0.992,0.553,0.235)
    RGB{Float64}(0.992,0.745,0.522)
    ]
    if length(sols) == 1
        labels = ["\$T_{f}\$"*labsuffix]
    else
        labels = ["\$T_{f$i}\$"*labsuffix for i in 1:length(sols)]
    end

    
    for (i, sol) in enumerate(sols)
        t_nd = range(0, sol.t[end-2], length=101)
        time = t_nd*u"hr"
        T = sol.(t_nd, idxs=2)*u"K"
        @series begin
            seriestype := :samplemarkers
            step := 20
            markershape --> :auto
            markersize --> 7
            seriescolor --> pal[i]
            label --> labels[i]
            return time, T
            end
        end
    end


    @doc raw"""
    modrftplot(sol, labsuffix=", model")
    modrftplot!(sol, labsuffix=", model")

Plot recipe for one solution to the lumped capacitance model.
This adds two series to the plot, with labels defaulting to `["T_f" "T_{vw}"] .* labsuffix`.
Since this is a recipe, any Plots.jl keyword arguments can be passed to modify the plot.
"""
modrftplot
@doc (@doc modrftplot) modrftplot!

@userplot ModRFTPlot
@recipe function f(tpmr::ModRFTPlot; labsuffix = ", model", evensample=true)
    sol = tpmr.args[1]
    if evensample
        t_nd = range(0, sol.t[end-2], length=31)
        step = 6
    else
        t_nd = sol.t[1:end-2]
        step = length(sol.t) ÷ 5
    end
    time = t_nd*u"hr"
    # color = palette(:Oranges_3)[end]
    color = RGB{Float64}(0.902,0.333,0.051)
    # Frozen temperature: tends to have a crazy time point at end
    @series begin
        seriestype := :samplemarkers
        step := step
        T = sol.(t_nd, idxs=2)*u"K"
        markershape --> :dtriangle
        markersize --> 7
        seriescolor --> color
        label --> "\$T_{f}\$"*labsuffix # Default label
        time, T
    end
    @series begin
        seriestype := :samplemarkers
        step := step
        T = sol.(t_nd, idxs=3)*u"K"
        markershape --> :utriangle
        markersize --> 7
        seriescolor --> color
        label --> "\$T_{vw}\$"*labsuffix # Default label
        linestyle := :dash
        time, T
    end
end

@doc raw"""
    tendplot(t_end)
    tendplot!(t_end)

Plot recipe that adds a labeled vertical line to the plot at time `t_end`. 
A default label and styling are applied, but these can be modified by keyword arguments as usual
"""
tendplot
@doc (@doc tendplot) tendplot!

@userplot tendPlot
@recipe function f(tp::tendPlot)
    t_end = tp.args[1]
    @series begin
        seriestype := :vline
        label --> "\$t_{end}\$"
        seriescolor --> :gray
        linestyle --> :dot
        linewidth --> 3
        return [ustrip(u"hr", t_end)]
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
        return ExpTfPlot((pdf.t_Tf, pdf.Tfs...))
    end
    if !ismissing(pdf.t_Tvw)
        @series begin
            return ExpTvwPlot((pdf.t_Tvw, pdf.Tvws...))
        end
    elseif !ismissing(pdf.Tvws)
        @series begin
            seriestype := :scatter
            return [pdf.t_Tf[end]], [pdf.t_Tvws]
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
            linewidth --> 0
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
    ys_cum = ys_cum ./ ys_cum[:,end]
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
        lumped_cap_rf_LC1(sol(ti), RF_params, ti)[2]
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
