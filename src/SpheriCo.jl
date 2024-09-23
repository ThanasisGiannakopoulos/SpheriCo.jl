module SpheriCo

using Parameters
using HDF5
using Interpolations

"""
Export function to build the spatial grid and structure for the
parameters used to build this grid.  These parameters are tunable, and
so the code can also run with an infalling outer boundary.
"""
export Grid, System
include("grid.jl")

"""
export functions defined in the following "include"
and called in classical and quantum modules
"""

export Dr_FD2, Drr_FD2, Dr_SBP2, KO_FD2, low_pass, project_on_smaller_grid
include("operators.jl")

"""
The time integrators are not exported here, but included at the level
of each submodule, explicitly and indepdendently.

TODO: find a way to export the time integrators at the level of the
main module, and just call them in the submodules, like e.g. the
operators above.
"""

# populate ghost points on the radial grid; needed for spatial derivatives
export even_ghosts, odd_ghosts
include("ghosts.jl")

# write various data
export write_data, write_quantum, write_data_r0, write_quantum_r0
export write_bilinears, write_correlators, write_debug, write_checkpoint
include("write_output.jl")

# useful function during evolution and in post-processing
export find_AH, energy, list_h5_files, Ricci_scalar
include("utils.jl")

include("classical/classical.jl")
include("quantum/quantum.jl")

end # module SpheriCo
