using Test
using SpheriCo

# Include other test files
include("test_ghosts.jl") 
include("test_FD2.jl")
include("test_SBP.jl")
include("test_time_integrators.jl")
include("Minkowski/test_ID_classic.jl")
include("Minkowski/test_constraints.jl")
include("Minkowski/test_rhs_classic.jl")
include("Minkowski/test_ID_quantum.jl")
include("Minkowski/test_rhs_quantum.jl")
include("Minkowski/test_rhs_quantum_trig_ID.jl")