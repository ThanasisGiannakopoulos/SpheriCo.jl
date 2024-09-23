# mutable structure with parameters to build the radial grid
@with_kw mutable struct Grid
    # discretization parameters
    Nr    :: Int
    r_max :: Float64
end

# the following structure saves the grid r and its spacing hr
struct System{T}
    r  :: Vector{T}
    hr :: T
end

""" 
the function System reads the Grid parameters, builds the grid, and
saves it in the structure System
"""
function System(grid::Grid)
    # -1 is for the last point on the right and -2 for ghosts on the left
    hr = grid.r_max / (grid.Nr-3)
    r  = range(-2*hr, grid.r_max, length=grid.Nr) 
    System{typeof(hr)}(r, hr)
end
