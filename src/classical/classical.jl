module classical

# export the function that runs the model, it is a function of Parameters
# i.e run_classical(p), with p=Param (see below for the Param structure)
export run_classical, Param

using SpheriCo

using Parameters
using Printf
using HDF5
using Dates
using Interpolations
using DelimitedFiles

import Base.Threads.@threads

include("parameters.jl")
include("rhs.jl")
include("../time_integrators.jl")
export constraints
include("constraints.jl")
include("ID.jl")
include("utils.jl")
include("evol.jl")

end #module
