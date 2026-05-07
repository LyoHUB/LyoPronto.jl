# Estimating Equipment Capability

The EC-CURT model, published by [Kazarin et al in 2021](https://link.springer.com/10.1208/s12249-021-02167-8), provides a compact way to get a first estimate of the equipment capability line. This equipment capability limit determines one edge of the design space for the primary drying stage of lyophilization.
This model is constructed by measuring equipment capability for 3 lyophilizers of various sizes, validating CFD calculations with those limits, then doing CFD calculations on a 3x3x3x2 experimental design: to evaluate on other lyophilizer geometries, a multilinear interpolation is used across this range of conditions. For full details, consult the article.

## Close translation

A numerical implementation in MATLAB for the EC-CURT model was provided in the supplementary information. The CFD data (chamber pressures at all conditions) were translated from there into Julia and provided as part of LyoPronto.jl.

```@docs
ECCURT.eq_cap_line
ECCURT.eq_cap_pressure
```

The original implementation uses MATLAB `interpn` with the `spline` method, which does a quadratic interpolation on dimensions with 3 sample points and linear interpolation on dimensions with 2 sample points. This Julia implementation uses Interpolations.jl, which only allows linear interpolations, so it will not return numerically identical results to the original MATLAB implementation.

## Newer application of same data

The original MATLAB implementation constructed lines by evaluating slope and intercept from the two higher mass flow rates at each geometry, neglecting the other two sampled points, and interpolating the slope and intercept values across lyophilizer geometries. To make full use of the reported data, we add a new implementation which interpolates the _pressures_ instead at all four mass flow rates, then regresses a new slope and intercept for the mass flow as a function of the interpolated pressures. 

This is provided here:

```@docs
ECCURT.eq_cap_line_new
```