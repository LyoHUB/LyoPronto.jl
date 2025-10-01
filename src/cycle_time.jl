# This code maybe belongs somewhere else but not sure where yet.
"""
    $(SIGNATURES)
    
Identify the end of primary drying using second derivative of Pirani pressure.


# Arguments
- `t`: time points
- `pch_pir`: Pirani pressure at given times
- `kind`: `Val(:der2)` for maximum of second derivative, `Val(:onoff)` for onset and offset.
- `window_width`: Width of the Savitzky-Golay filter window (default: 91)
- `tmin`: Minimum time for analysis (default: 0 hours, unitful)
- `tmax`: Maximum time for analysis (default: Inf hours, unitful)

If a single argument `data` is passed instead of `t` and `pch_pir`, the time and Pirani pressure are assumed to be 
accessible as either `data.t` and `data.pch_pir` or `data[:t]` and `data[:pch_pir]`.

# Approach

Use a Savitzky-Golay filter to smooth and take derivatives.
For the onset-offset, find linear intersects of tangent at inflection point with maximum and 
minimum values of pressure. Those intersects are taken as onset and offset.
For the second derivative, identify the maximum of the second derivative.
"""
function identify_pd_end(t, pch_pir, kind; window_width=91, tmin=0u"hr", tmax=Inf*u"hr")
    if kind isa Val{:der2}
        pch_pir_der2 = savitzky_golay(pch_pir, window_width, 3, deriv=2).y
        t_end = t[argmax(pch_pir_der2[tmin .< t .< tmax])] + tmin
        return t_end
    elseif kind isa Val{:onoff}
        pch_pir_sm = savitzky_golay(pch_pir, window_width, 3, deriv=0).y
        dt = sum(diff(t)) / (length(t)-1)
        pch_pir_der1 = savitzky_golay(pch_pir, window_width, 3, deriv=1, rate=1/dt).y
        ti_inds = tmin .< t .< tmax
        i_mid = argmin(pch_pir_der1[ti_inds]) + searchsortedfirst(t, tmin)
        t_mid = t[i_mid]
        dp_mid = pch_pir_der1[i_mid]
        p_mid = pch_pir_sm[i_mid]
        p_max = maximum(pch_pir_sm[ti_inds])
        p_min = minimum(pch_pir_sm[ti_inds])
        t_onset = (p_max - p_mid)/dp_mid + t_mid
        t_offset = (p_min - p_mid)/dp_mid + t_mid
        return (t_onset, t_offset)
    else
        error("Unknown kind $kind. Use :der2 or :onoff")
    end
end
function identify_pd_end(data, kind; window_width=91, tmin=0u"hr", tmax=Inf*u"hr")
    if hasproperty(data, :t) && hasproperty(data, :pch_pir)
        t = data.t
        pch_pir = data.pch_pir
    else
        t = data[:t]
        pch_pir = data[:pch_pir]
    end
    return identify_pd_end(t, pch_pir, kind; window_width, tmin, tmax)
end
