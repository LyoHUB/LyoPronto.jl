export get_vial_radii, get_vial_mass, get_vial_shape, make_outlines
# Get dimensions from here: https://www.schott.com/en-us/products/vials/-/media/project/onex/products/v/vials/application-variants/schott-brochure-schott-vials-english-20092017.pdf

# Added to "vial_sizes.csv", read in here
const VIAL_DIMS = CSV.File(srcdir("vial_sizes.csv"))

"""
    get_vial_radii(vialsize::String)

Return inner and outer radius for passed ISO vial size.

Uses a table provided by Schott, stored internally in a CSV.
"""
function get_vial_radii(vialsize::String)
    alldims = [row for row in VIAL_DIMS if row.Size == vialsize]
    if length(alldims) != 1
        @error "bad vial size passed" vialsize
    end
    alldims = alldims[1] # Extract object corresponding to row of table
    rad_o = alldims.d1 / 2 * u"mm"
    rad_i = rad_o - alldims.s1 * u"mm"
    return rad_i, rad_o
end

"""
    get_vial_thickness(vialsize::String)

Return vial wall thickness for given ISO vial size.

Uses a table provided by Schott, stored internally in a CSV.
"""
function get_vial_thickness(vialsize::String)
    alldims = [row for row in VIAL_DIMS if row.Size == vialsize]
    if length(alldims) != 1
        @error "bad vial size passed" vialsize
    end
    alldims = alldims[1] # Extract object corresponding to row of table
    thickness = alldims.s1*u"mm"
    return thickness
end

function get_vial_mass(vialsize::String)
    alldims = [row for row in VIAL_DIMS if row.Size == vialsize]
    if length(alldims) != 1
        @error "bad vial size passed" vialsize
    end
    alldims = alldims[1] # Extract object corresponding to row of table
    thickness = alldims.mass*u"g"
    return thickness
end

function get_vial_shape(vialsize::String)
    alldims = [row for row in VIAL_DIMS if row.Size == vialsize]
    if length(alldims) != 1
        @error "bad vial size passed" vialsize
    end
    alldims = alldims[1] # Extract object corresponding to row of table
    rad_o = alldims.d1 / 2 * u"mm"
    rad_i = rad_o - alldims.s1 * u"mm"
    bot_thick = alldims.s2*u"mm"
    full_height = alldims.h1*u"mm"
    curve_height = full_height - alldims.h3*u"mm"
    barrel_height = alldims.h2*u"mm"
    neck_inner = alldims.d4/2*u"mm"
    neck_outer = alldims.d3/2*u"mm"
    neck_curve = alldims.r1*u"mm"
    return @dict rad_i rad_o bot_thick barrel_height curve_height full_height neck_inner neck_outer neck_curve
end


"""
    make_outlines(dims, Vfill)

Return `Plots.Shape`s for vial and fill volume, with Unitful dimensions, for given vial dimensions.
Convenience function for making some figures.
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
    return Shape(vpoints), Shape(fpoints)
end
