module quantum

# copied from Jecco.jl/src/AdS5_3_1/AdS5_3_1.jl
# needed to set OPENBLAS_NUM_THREADS=1 automatically
using LinearAlgebra

# always set the number of BLAS threads to 1 upon loading the module. by default
# it uses a bunch of them and we don't want that since they trample over each
# other when solving the nested systems equations. it's much better to thread
# over the loop. see also the discussion here:
# https://github.com/JuliaLang/julia/issues/33409
#
# this saves us the tedious task of always setting OMP_NUM_THREADS=1 before
# launching julia.
function __init__()
    LinearAlgebra.BLAS.set_num_threads(1)
    nothing
end

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
