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
    tplotexperimental(time, T1, [T2, ...])
    tplotexperimental!(time, T1, [T2, ...])

Plot recipe for one or more experimentally measured product temperatures, all at same times.
This recipe adds one series for each passed temperature series, so pass labels as appropriate.
"""
tplotexperimental
@doc (@doc tplotexperimental) tplotexperimental!

@userplot TPlotExperimental
@recipe function f(tpe::TPlotExperimental)
    time, Ts... = tpe.args
    step = size(time, 1) ÷ 10
    n = size(Ts, 1)
    # pal = palette(:Blues_5,)[end:-1:begin] # Requires Plots for palette, so hard-code this default
    pal = [RGB{Float64}(0.031,0.318,0.612), 
           RGB{Float64}(0.192,0.51,0.741), 
           RGB{Float64}(0.42,0.682,0.839), 
           RGB{Float64}(0.741,0.843,0.906)]
    
    for (i, T) in enumerate(Ts)
        @series begin
            seriestype := :samplemarkers
            step := step
            offset := step÷n *(i-1) + 1
            # markershape --> :auto
            markersize --> 7
            # label := labels[i]
            minlen = min(length(time), length(T))
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
    tplotexpvw(time, temperature; trim)
    tplotexpvw!(time, temperature; trim)

Plot recipe for a set of experimentally measured vial wall temperatures.
This recipe adds only one series to the plot.
`trim` is an integer, indicating how many points to skip at a time, so that 
the dotted line looks dotted even with noisy data.
"""
tplotexpvw
@doc (@doc tplotexpvw) tplotexpvw!

@userplot TPlotExpVW
@recipe function f(tpev::TPlotExpVW; trim=1)
    time, T = tpev.args
    time_trim = time[begin:trim:end]
    T_trim = T[begin:trim:end]
    step = size(time_trim, 1) ÷ (10)
    # color = palette(:Blues_5)[end]
    color = RGB{Float64}(0.031,0.318,0.612)
    
    @series begin
        seriestype := :samplemarkers
        step := step
        markershape --> :rect
        markersize --> 7
        seriescolor --> color
        linestyle := :dash
        # x := time
        # y := T
        time_trim, T_trim
    end
end

@doc raw"""
    tplotmodelconv(sols)
    tplotmodelconv!(sols)

Plot recipe for one or multiple solutions to the Pikal model, e.g. the output of [`gen_sol_conv_dim`](@ref LyoPronto.gen_sol_conv_dim).
This adds one series to the plot for each passed solution, so pass as many labels (e.g. `["Tf1" "Tf2"]`) to this plot call as solutions to add labels to the legend.
"""
tplotmodelconv
@doc (@doc tplotmodelconv) tplotmodelconv!


@userplot TPlotModelConv
@recipe function f(tpmc::TPlotModelConv)
    sols = tpmc.args
    # pal = palette(:Oranges_4).colors[end:-1:begin+1] # Requires Plots as dependency...
    pal = [
    RGB{Float64}(0.851,0.278,0.004)
    RGB{Float64}(0.992,0.553,0.235)
    RGB{Float64}(0.992,0.745,0.522)
    ]
    
    for (i, sol) in enumerate(sols)
        @series begin
            t_nd = range(0, sol.t[end-2], length=101)
            time = t_nd*u"hr"
            T = sol.(t_nd, idxs=2)*u"K"
            seriestype := :samplemarkers
            step := 20
            # offset := step÷n *(i-1) + 1
            markershape --> :auto
            markersize --> 7
            seriescolor --> pal[i]
            # x := time
            # y := T
            time, T
        end
    end
end


@doc raw"""
    tplotmodelrf(sol)
    tplotmodelrf!(sol)

Plot recipe for one solution to the lumped capacitance model, e.g. the output of [`gen_sol_rf_dim`](@ref LyoPronto.gen_sol_rf_dim).
This adds two series to the plot, so pass two labels (e.g. `["Tf" "Tvw"]`) to this plot call to add labels to the legend.
"""
tplotmodelrf
@doc (@doc tplotmodelrf) tplotmodelrf!

@userplot TPlotModelRF
@recipe function f(tpmr::TPlotModelRF)
    sol = tpmr.args[1]
    t_nd = range(0, sol.t[end-2], length=31)
    time = t_nd*u"hr"
    # color = palette(:Oranges_3)[end]
    color = RGB{Float64}(0.902,0.333,0.051)
    # Frozen temperature: tends to have a crazy time point at end
    @series begin
        # time = sol.t[1:end-2]*u"hr"
        seriestype := :samplemarkers
        step := 6
        T = sol.(t_nd, idxs=2)*u"K"
        markershape --> :dtriangle
        markersize --> 7
        seriescolor --> color
        time, T
    end
    @series begin
        # time = sol.t*u"hr"
        # T = sol[3,:]*u"K"
        seriestype := :samplemarkers
        step := 6
        T = sol.(t_nd, idxs=3)*u"K"
        markershape --> :utriangle
        markersize --> 7
        seriescolor --> color
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
        label --> "\$t_\\text{end}\$"
        color --> :gray
        ls --> :dot
        return [ustrip(u"hr", t_end)]
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
    # @info "check" sol integ_values

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
    return TPlotExperimental((pdf.t_Tf, pdf.Tfs...))
end
@recipe function f(pdf::PrimaryDryFit, vw::Val{true})
    if ismissing(pdf.t_Tvw)
        return [pdf.t_Tf[end]], [pdf.t_Tvws]
    else
        return TPlotExpVW((pdf.t_Tvw, pdf.Tvws[1]))
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


# function qplotrf(sol, RF_params; kw...)
#     Qcontrib = map(sol.t) do ti
#         lumped_cap_rf(sol(ti), RF_params, ti, energy_output=true)[2]
#     end
#     Qcontrib = hcat(Qcontrib...)
#     Qsub = Qcontrib[1,:]
#     Qshf = Qcontrib[2,:]
#     Qvwf = Qcontrib[3,:]
#     QRFf = Qcontrib[4,:]
#     QRFvw = Qcontrib[5,:]
#     # names = ["sub", "sh-f", "vw-f", "RF-f", "RF-vw", "sh-vw"]
#     names = ["RF-f", "vw-f", "sh-f"]
#     labs = ["\$Q_\\text{$nm}\$" for nm in names]
#     pl = plot(u"hr", u"W", kw...)
#     areastackplot!(sol.t, QRFf, Qvwf, Qshf, labels=permutedims(labs))
#     plot!(sol.t, Qshf.+Qvwf.+QRFf, c=:black, label="Total")
#     plot!(xlabel="Time [hr]", ylabel="Heating [W]")
#     return pl
# end

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
