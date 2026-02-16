# Get dimensions from here: https://www.schott.com/en-us/products/vials/-/media/project/onex/products/v/vials/application-variants/schott-brochure-schott-vials-english-20092017.pdf

# That table is compact enough to include here as a commented-out CSV

# Size,overflow,overflow_tol,a,d1,d1_tol,d2,d3,d4,h1,h1_tol,h2,h3,h3_tol,r1,r2,s1,s1_tol,s2,t,mass
# 2R,4,0.5,1,16,0.15,13,10.5,7,35,0.5,22,8,0.5,2.5,1.5,1,0.04,0.6,0.7,4.4
# 4R,6,0.5,1,16,0.15,13,10.5,7,45,0.5,32,8,0.5,2.5,1.5,1,0.04,0.6,0.7,5.7
# 6R,10,0.5,1.2,22,0.2,20,16.5,12.6,40,0.5,26,8.5,0.5,3.5,2,1,0.04,0.7,0.7,7.9
# 8R,11.5,0.5,1.2,22,0.2,20,16.5,12.6,45,0.5,31,8.5,0.5,3.5,2,1,0.04,0.7,0.7,8.7
# 10R,13.5,1,1.2,24,0.2,20,16.5,12.6,45,0.5,30,9,0.5,4,2,1,0.04,0.7,0.7,9.5
# 15R,19,1,1.2,24,0.2,20,16.5,12.6,60,0.5,45,9,0.5,4,2,1,0.04,0.7,0.7,12
# 20R,26,1.5,1.5,30,0.25,20,17.5,12.6,55,0.7,35,10,0.75,5.5,2.5,1.2,0.05,0.7,1,16.2
# 25R,32.5,1.5,1.5,30,0.25,20,17.5,12.6,65,0.7,45,10,0.75,5.5,2.5,1.2,0.05,0.7,1,18.9
# 30R,37.5,1.5,1.5,30,0.25,20,17.5,12.6,75,0.7,55,10,0.75,5.5,2.5,1.2,0.05,0.7,1,21.9
# 50R,62,4,2.5,40,0.4,20,17.5,12.6,73,0.75,49,10,0.75,6,4,1.5,0.07,0.9,1.5,34.5
# 100R,123,7,3.5,47,0.5,20,17.5,12.6,100,0.75,75,10,0.75,6.5,4,1.7,0.07,0.9,1.5,60

# This is transformed mostly-programmatically into the following DictTable:

VIAL_DIMS = DictTable(
    Size = ["2R", "4R", "6R", "8R", "10R", "15R", "20R", "25R", "30R", "50R", "100R"],
    overflow = [4.0, 6.0, 10.0, 11.5, 13.5, 19.0, 26.0, 32.5, 37.5, 62.0, 123.0]*u"mL",
    overflow_tol = [0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.5, 1.5, 1.5, 4.0, 7.0]*u"mL",
    a = [1.0, 1.0, 1.2, 1.2, 1.2, 1.2, 1.5, 1.5, 1.5, 2.5, 3.5]*u"mm",
    d1 = [16, 16, 22, 22, 24, 24, 30, 30, 30, 40, 47]*u"mm",
    d1_tol = [0.15, 0.15, 0.2, 0.2, 0.2, 0.2, 0.25, 0.25, 0.25, 0.4, 0.5]*u"mm",
    d2 = [13, 13, 20, 20, 20, 20, 20, 20, 20, 20, 20]*u"mm",
	d3 = [10.5, 10.5, 16.5, 16.5, 16.5, 16.5, 17.5, 17.5, 17.5, 17.5, 17.5]*u"mm",
	d4 = [7.0, 7.0, 12.6, 12.6, 12.6, 12.6, 12.6, 12.6, 12.6, 12.6, 12.6]*u"mm",
	h1 = [35, 45, 40, 45, 45, 60, 55, 65, 75, 73, 100]*u"mm",
	h1_tol = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.7, 0.7, 0.7, 0.75, 0.75]*u"mm",
	h2 = [22, 32, 26, 31, 30, 45, 35, 45, 55, 49, 75]*u"mm",
	h3 = [8.0, 8.0, 8.5, 8.5, 9.0, 9.0, 10.0, 10.0, 10.0, 10.0, 10.0]*u"mm",
	h3_tol = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75]*u"mm",
	r1 = [2.5, 2.5, 3.5, 3.5, 4.0, 4.0, 5.5, 5.5, 5.5, 6.0, 6.5]*u"mm",
	r2 = [1.5, 1.5, 2.0, 2.0, 2.0, 2.0, 2.5, 2.5, 2.5, 4.0, 4.0]*u"mm",
	s1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.2, 1.2, 1.2, 1.5, 1.7]*u"mm",
	s1_tol = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.07, 0.07]*u"mm",
	s2 = [0.6, 0.6, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9]*u"mm",
	t = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 1.0, 1.0, 1.0, 1.5, 1.5]*u"mm",
	mass = [4.4, 5.7, 7.9, 8.7, 9.5, 12.0, 16.2, 18.9, 21.9, 34.5, 60.0]*u"g",
)

const VIAL_DIMS_SOURCE_DOC = """
Uses a table from a SCHOTT manual, stored internally in a CSV.
"""

function select_size(vialsize::String)
    if vialsize ∉ VIAL_DIMS.Size
        error("bad vial size passed: $vialsize")
    end
    return VIAL_DIMS[vialsize]
end

"""
    $(SIGNATURES)

Return inner and outer radius for passed ISO vial size.

$(VIAL_DIMS_SOURCE_DOC)
"""
function get_vial_radii(vialsize::String)
    alldims = select_size(vialsize)
    rad_o = alldims.d1 / 2
    rad_i = rad_o - alldims.s1
    return rad_i, rad_o
end

"""
    $(SIGNATURES)

Return vial wall thickness for given ISO vial size.

$(VIAL_DIMS_SOURCE_DOC)
"""
get_vial_thickness(vialsize::String) = select_size(vialsize).s1

"""
    $(SIGNATURES)

Return vial mass for given ISO vial size.

$(VIAL_DIMS_SOURCE_DOC)
"""
get_vial_mass(vialsize::String) = select_size(vialsize).mass

"""
    get_vial_shape(vialsize::String)

Return a NamedTuple with a slew of vial dimensions, useful for drawing the shape of the vial with [`make_outlines`](@ref).

$(VIAL_DIMS_SOURCE_DOC)
"""
function get_vial_shape(vialsize::String)
    alldims = select_size(vialsize)
    rad_o = alldims.d1 / 2
    rad_i = rad_o - alldims.s1
    bot_thick = alldims.s2
    full_height = alldims.h1
    curve_height = full_height - alldims.h3
    barrel_height = alldims.h2
    neck_inner = alldims.d4/2
    neck_outer = alldims.d3/2
    neck_curve = alldims.r1
    dims = (; rad_i, rad_o, bot_thick, barrel_height, curve_height, full_height, neck_inner, neck_outer, neck_curve)
    return dims
end


@doc raw"""
    make_outlines(dims, Vfill)

Return a sequence of points (ready to be made into `Plots.Shape`s for the vial and fill volume, with Unitful dimensions, for given vial dimensions.

This is a convenience function for making figures illustrating fill depth.
"""
function make_outlines(dims, Vfill)
    (; rad_o, rad_i, bot_thick, neck_inner, neck_outer, curve_height, barrel_height, full_height) = dims

    vpoints = [ 
        (-rad_o, 0.0*u"mm"),
        ( rad_o, 0.0*u"mm"),
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
        (-rad_o, 0.0*u"mm"),
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
