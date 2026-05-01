module ECCURT
using Accessors
using Unitful: ustrip, @u_str, Quantity
using DocStringExtensions
using Interpolations: interpolate, Gridded, Linear, extrapolate, Line

public eq_cap_pressure, eq_cap_line

const M_DOT = [0.1291, 0.4644, 0.7776, 1.1772]

const D_SAMPLE = [59.6, 98.0, 304.8] # mm
const DA_SAMPLE = [2.24, 6.49, 21.75] # diameter / valve thickness (unitless)
const L_SAMPLE = reverse([890.0, 441.0, 100.0]) # Reversed order, so increasing
const VOLUME_SAMPLE = [0.092, 0.44]

# First dimension: mdot sample points, written vectors in original MatLab
# Second dimension: chamber size, f & g in original MatLab
# Third dimension: 9 points for each spool (3 valve thickness?, 3 length?)
# Fourth dimension: spool size (S20, LS3, FDX) in original MatLab
const PCH = [
    # # LS3 chamber left half, S20 chamber right half
    # S20 spool
    227.9; 557.5; 865.4; 1258.2;; 223.8; 556.8; 867.8; 1264.6;;;
    110.7; 261.1; 400.0; 579.2;; 110.3; 260.2; 400.2; 578.8;;;
    13.4; 32.5; 48.4; 68.0;; 12.7; 30.9; 46.15; 65.0;;;;
    175.2; 374.3; 560.2; 797.5;; 177.9; 390.8; 589.6; 843.3;;;
    79.0; 166.6; 248.3; 352.6;; 64.9; 159.3; 247.5; 360.0;;;
    11.9; 22.8; 33.0; 46.0;; 11.7; 21.1; 30.0; 41.3;;;;
    169.4; 351.7; 522.0; 740.0;; 173.2; 358.1; 530.9; 751.3;;;
    75.0; 154.3; 228.0; 322.0;; 75.9; 154.6; 228.0; 321.7;;;
    8.9; 20.9; 30.6; 42.5;; 9.8; 18.69; 27.0; 37.6;;;;;
    # LS3 spool
    241.3; 635.8; 1004.1; 1474.2;; 248.5; 644.4; 1014.1; 1485.8;;;
    93.7; 250.5; 390.4; 569.3;; 93.7; 250.6; 390.9; 569.8;;;
    15.1; 31.8; 47.4; 67.3;; 14.4; 30.4; 45.4; 64.5;;;;
    158.4; 377.2; 581.6; 842.4;; 168.4; 399.1; 614.7; 889.7;;;
    67.0; 151.7; 230.9; 331.9;; 66.5; 151.1; 230.1; 331.0;;;
    10.9; 21.5; 31.5; 44.2;; 14.4; 22.5; 30.0; 39.6;;;;
    148.6; 344.3; 527.1; 760.3;; 151.7; 351.1; 537.4; 775.1;;;
    60.9; 136.7; 207.6; 298.0;; 60.5; 136.1; 206.6; 296.6;;;
    10.0; 19.8; 29.0; 40.7;; 8.28; 16.93; 25.0; 35.3;;;;;
    # FDX spool
    203.8; 605.4; 980.4; 1458.9;; 158.5; 571.2; 956.7; 1448.6;;;
    88.0; 239.5; 381.0; 561.5;; 88.5; 239.4; 380.4; 560.2;;;
    14.3; 31.2; 47.0; 67.2;; 13.3; 29.5; 44.7; 64.0;;;;
    116.4; 331.7; 532.9; 789.5;; 118.6; 338.2; 543.2; 804.8;;;
    50.2; 133.8; 212.0; 311.6;; 50.5; 133.4; 210.9; 309.7;;;
    9.9; 20.6; 30.6; 43.4;; 9.3; 18.6; 27.3; 38.5;;;;
    107.2; 295.1; 470.6; 694.5;; 109.1; 299.9; 478.1; 705.5;;;
    43.6; 116.8; 185.2; 272.4;; 43.9; 116.4; 184.1; 270.5;;;
    9.2; 18.9; 28.0; 39.5;; 7.3; 15.8; 23.8; 33.9;;;;;
]

struct ECLine{K,B}
    k::K
    b::B
end
function (line::ECLine)(p)
    return line.k * p + line.b
end

# Original MATLAB code: 
# for i=1:27
#     A(i).k = (m_dot(3)-m_dot(4))/(Pch(i).f(3)-Pch(i).f(4));
#     A(i).b = m_dot(3)-A(i).k*Pch(i).f(3);
# end
# Dimensions here: d, da, l?
const Ak = [(M_DOT[3] - M_DOT[4])/(PCH[3,1,i,j,k] - PCH[4,1,i,j,k]) for i in axes(PCH,3), j in axes(PCH,4), k in axes(PCH, 5)]
const Ab = [M_DOT[3] - Ak[i,j,k]*PCH[3,1,i,j,k] for i in axes(PCH,3), j in axes(PCH,4), k in axes(PCH, 5)]

# for i=1:27
#     A(i+27).k = (m_dot(3)-m_dot(4))/(Pch(i).g(3)-Pch(i).g(4));
#     A(i+27).b = m_dot(3)-A(i+27).k*Pch(i).g(3);
# end
# Dimensions here: d, da, l?
const Bk = [(M_DOT[3] - M_DOT[4])/(PCH[3,2,i,j,k] - PCH[4,2,i,j,k]) for i in axes(PCH,3), j in axes(PCH,4), k in axes(PCH, 5)]
const Bb = [M_DOT[3] - Bk[i,j,k]*PCH[3,2,i,j,k] for i in axes(PCH,3), j in axes(PCH,4), k in axes(PCH, 5)]

# Last dimension from stacking: v
# All dimensions: d, da, l, v
# Last two dimensions are inverted from MATLAB order
# l dimension reversed to put sample points in increasing order
const alpha = stack([Ak[:,:,end:-1:begin], Bk[:,:,end:-1:begin]])
const beta = stack([Ab[:,:,end:-1:begin], Bb[:,:,end:-1:begin]])

const alpha_interp = interpolate((D_SAMPLE, DA_SAMPLE, L_SAMPLE, VOLUME_SAMPLE), alpha, Gridded(Linear()))
const alpha_ext = extrapolate(alpha_interp, Line())
const beta_interp = interpolate((D_SAMPLE, DA_SAMPLE, L_SAMPLE, VOLUME_SAMPLE), beta, Gridded(Linear()))
const beta_ext = extrapolate(beta_interp, Line())

"""
    $(SIGNATURES)

Compute the equipment capability line for given geometry parameters.
The resulting object can be called with a pressure in mTorr to get the corresponding mass flow rate in kg/hr.
The line's slope and intercept can be accessed with `line.k` and `line.b`, respectively.

If provided with `Unitful.Quantity` inputs, the geometry parameters will be converted to the
expected units (mm for dimensions, m^3 for volume) before calculation, and the result will 
have Unitful units as well. If provided with plain floats, plain floats will be returned.

Dimensions:
- `d`: Diameter of the duct (spool) from chamber to condenser [mm]
- `vt`: Thickness of the butterfly valve in the duct (spool) from chamber to condenser [mm]
- `l`: Length of the duct (spool) [mm]
- `volume`: Volume of the product chamber [m^3]
"""
function eq_cap_line(d, vt, l, volume)
    if ~(D_SAMPLE[1] <= d <= D_SAMPLE[end]) || ~(DA_SAMPLE[1] <= d/vt <= DA_SAMPLE[end]) || ~(L_SAMPLE[1] <= l <= L_SAMPLE[end]) || ~(VOLUME_SAMPLE[1] <= volume <= VOLUME_SAMPLE[end])
        @warn "Input geometry parameters are outside the range of the original data, so extrapolation is being used. Results may be inaccurate."
    end
    alpha = alpha_ext(d, d/vt, l, volume)
    beta = beta_ext(d, d/vt, l, volume)

    return ECLine(alpha, beta)
end
function eq_cap_line(d::Quantity, vt::Quantity, l::Quantity, volume::Quantity)
    ndim = eq_cap_line(ustrip(u"mm", d), ustrip(u"mm", vt), 
        ustrip(u"mm", l), ustrip(u"m^3", volume))
    @reset ndim.k *= u"kg/hr/mTorr"
    @reset ndim.b *= u"kg/hr"
    return ndim
end

"""
    $(SIGNATURES)

Call `eq_cap_line` to get the equipment capability line for the given geometry parameters, 
then evaluate that line with mass flow rate `m` in kg/hr to get minimum controllable pressure in mTorr.
"""
function eq_cap_pressure(m, D, valve_thickness, L, volume)
    line = eq_cap_line(D, valve_thickness, L, volume)
    return (m - line.b) / line.k
end


end # module ECCURT