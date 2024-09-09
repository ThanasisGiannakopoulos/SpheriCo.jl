module SpheriCo

# export structure for tunable parameters
export Grid, System

using Parameters
using HDF5
using Interpolations

# to construct the radial grid; needed for derivative operators
export System
include("grid.jl")

"""
export functions defined in the following "include
and called in classical and quantum modules
"""

export Dr_FD2, Drr_FD2, Dr_SBP2, KO_FD2, low_pass, project_on_smaller_grid
include("operators.jl")

"""
TODO export the rhs function F! such that I can export the time integrators here
 and feed them the different F! later in the submodules
"""

# """
# empty function F!(t::Float64, vf::Array, v::Array, sys::System, p::Param)
# needed for the time integrators. It is specified later in the modules classical and quantum
# """
# function F!(t::Float64, vf::Array, v::Array, sys::System, p::Param) end
# export RK4, AB3!
# include("time_integrators.jl")

export even_ghosts, odd_ghosts
include("ghosts.jl")

export write_data, write_quantum, write_data_r0, write_quantum_r0
export write_bilinears, write_correlators, write_debug, write_checkpoint
include("write_output.jl")

export find_AH, energy, list_h5_files, Ricci_scalar #, constraints_sbp21
include("utils.jl")

include("classical/classical.jl")
include("quantum/quantum.jl")

end # module SpheriCo
