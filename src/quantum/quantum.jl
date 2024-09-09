module quantum

# export the function that runs the model, it is a function of Parameters
# i.e run_classical(p), with p=Param (see below for the Param structure)
export run_quantum, Param

using SpheriCo

using Parameters
using Printf
using HDF5
using Dates
using Interpolations
using DelimitedFiles
using SpecialFunctions

import Base.Threads.@threads

# parameters for the model
# the grid is given in SpheriCo modul as mutable structure
include("parameters.jl")

# from SpheriCo for time integration;
# TODO: need to find a proper way to extract the functions
include("../time_integrators.jl")

# for the classical part of the simulation
# Hamiltonian & Momentum constraints
include("../classical/constraints.jl")

# classical ID, real-valued
include("../classical/ID.jl")
# quantum ID, complex-valued
include("ID.jl")

# bilinears and correlators
include("bilinears.jl")
include("correlators.jl")

include("rhs.jl")
include("utils.jl")
include("evol.jl")

end #module
