# mutable structure with parameters to build the radial grid
@with_kw mutable struct Grid
    # discretization parameters
    Nr    :: Int
    r_max :: Float64
end

struct System{T}
    r     :: Vector{T}
    hr    :: T
end
function System(grid::Grid)
    # -1 is for the last point on the right and -2 for ghosts on the left
    hr = grid.r_max / (grid.Nr-3)
    r  = range(-2*hr, grid.r_max, length=grid.Nr) 
    System{typeof(hr)}(r, hr)
end
