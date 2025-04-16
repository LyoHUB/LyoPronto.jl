export get_vial_radii, get_vial_mass, get_vial_shape, make_outlines
# Get dimensions from here: https://www.schott.com/en-us/products/vials/-/media/project/onex/products/v/vials/application-variants/schott-brochure-schott-vials-english-20092017.pdf

# Added to "vial_sizes.csv", read in here
const VIAL_DIMS = CSV.File((@__DIR__) * raw"/vial_sizes.csv")

function select_size(vialsize::String)
    alldims = filter(x->x.Size == vialsize, VIAL_DIMS)
    if length(alldims) != 1
        @error "bad vial size passed" vialsize
    end
    alldims = alldims[1] # Extract object corresponding to row of table
end

@doc raw"""
    get_vial_radii(vialsize::String)

Return inner and outer radius for passed ISO vial size.

Uses a table provided by Schott, stored internally in a CSV.
"""
function get_vial_radii(vialsize::String)
    alldims = select_size(vialsize)
    rad_o = alldims.d1 / 2 * u"mm"
    rad_i = rad_o - alldims.s1 * u"mm"
    return rad_i, rad_o
end

@doc raw"""
    get_vial_thickness(vialsize::String)

Return vial wall thickness for given ISO vial size.

Uses a table provided by Schott, stored internally in a CSV.
"""
function get_vial_thickness(vialsize::String)
    alldims = select_size(vialsize)
    thickness = alldims.s1*u"mm"
    return thickness
end

@doc raw"""
    get_vial_mass(vialsize::String)

Return vial mass for given ISO vial size.

Uses a table provided by Schott, stored internally in a CSV.
"""
function get_vial_mass(vialsize::String)
    alldims = select_size(vialsize)
    mass = alldims.mass*u"g"
    return mass
end

@doc raw"""
    get_vial_shape(vialsize::String)

Return a Dict{Symbol, } with a slew of vial dimensions, useful for drawing the shape of the vial with [`make_outlines`](@ref).
"""
function get_vial_shape(vialsize::String)
    alldims = select_size(vialsize)
    rad_o = alldims.d1 / 2 * u"mm"
    rad_i = rad_o - alldims.s1 * u"mm"
    bot_thick = alldims.s2*u"mm"
    full_height = alldims.h1*u"mm"
    curve_height = full_height - alldims.h3*u"mm"
    barrel_height = alldims.h2*u"mm"
    neck_inner = alldims.d4/2*u"mm"
    neck_outer = alldims.d3/2*u"mm"
    neck_curve = alldims.r1*u"mm"
    dims = Dict{Symbol, typeof(rad_o)}()
    @pack! dims = rad_i, rad_o, bot_thick, barrel_height, curve_height, full_height, neck_inner, neck_outer, neck_curve
    return dims
end


@doc raw"""
    make_outlines(dims, Vfill)

Return a sequence of points (ready to be made into `Plots.Shape`s for the vial and fill volume, with Unitful dimensions, for given vial dimensions.

This is a convenience function for making figures illustrating fill depth.
"""
function make_outlines(dims, Vfill)
    @unpack rad_o, rad_i, bot_thick, neck_inner, neck_outer, curve_height, barrel_height, full_height = dims

    vpoints = [ 
        (-rad_o, 0*u"mm"),
        ( rad_o, 0*u"mm"),
        ( rad_o, barrel_height),
        ( neck_outer, curve_height),
        ( neck_outer, full_height),
        ( neck_inner, full_height),
        ( neck_inner, curve_height),
        ( rad_i, barrel_height),
        ( rad_i, bot_thick),
        (-rad_i, bot_thick),
        (-rad_i, barrel_height),
        (-neck_inner, curve_height),
        (-neck_inner, full_height),
        (-neck_outer, full_height),
        (-neck_outer, curve_height),
        (-rad_o, barrel_height),
        (-rad_o, 0*u"mm"),
    ]
    fheight = Vfill / (π*rad_i^2)
    fpoints = [
        (-rad_i, bot_thick),
        (rad_i, bot_thick),
        (rad_i, bot_thick+fheight),
        (-rad_i, bot_thick+fheight),
        (-rad_i, bot_thick),
    ]
    return vpoints, fpoints
end
