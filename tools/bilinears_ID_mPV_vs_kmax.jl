"""
With this script you can check the effect of quantum modes in the
backreaction, at the level of the initial data. It assumes Minkowski
space, and by tuning the number of quantum modes included (tune k, l
quantum numbers), you can see how close to Minkowski is the source
term in the equations of motion for the extrinsic curvature components
K and KB.

The goal is to check what is the minimum kmax for which the
backreaction is reasonable, for different values of mPV. By reasonable
here we mean that the matter content for the extrinsic curvature rhs
tends to zero with increasing number of kmax. We only focus on r=0 and
we want to explore the interplay with mPV. See
https://arxiv.org/pdf/2111.11400 for how the interplay between kmax
and lmax affects the backreaction.
"""

# load the needed packages and functions
using HDF5
using LaTeXStrings
using Plots ; pythonplot()
using DelimitedFiles
using SpecialFunctions
using FunctionZeros
using Interpolations

using SpheriCo
using SpheriCo.quantum
include("../src/classical/ID.jl")
import Base.Threads.@threads
include("../src/quantum/ID.jl")
include("../src/quantum/bilinears.jl")

# you can parallelize
println("Running with number of threads = ", Threads.nthreads())

out_dir = "./bilinears_ID/"

# Check if the directory exists, and create it if it doesn't
if !isdir(out_dir)
    mkdir(out_dir)
    println("Directory created: $out_dir")
else
    println("Directory already exists: $out_dir")
end

# function that calculates the bilinears, D is used to tune the number pf grid points, kmax, lmax are quantum numbers
# returns v_classic, v_quantum, bln, ρ_all, SA_all, SB_all, jA_all, rr
function calculate_bln(D, lmax, kmax, mPV)

    # convention for the physical Planck mass M_P^2, 1 or 1/(8pi)
    Mp2_phys = 1#/(8*π)
    # what the code sees in 1/M_P^2
    oMp2 = 1/(Mp2_phys - mPV^2 * log(2^4/3^3) / (12*4*π^2))
    
    # the theoretical value for the cosmological constant, when backreaction is included
    CC_theory = -log(3^9/2^16)*(mPV^4)/(8*(2*π)^2)
    
    # the denominator of dk says how far in r the backreaction is good (approx up until the denom of dk)
    dk = π/30
    
    # convention, quantum
    hbar = 1.0
    c = 1.0
    
    # radial grid
    g = Grid(
        # discretization parameters
        Nr         = (128)*2^D + 1 + 2,
        r_max      = 30.0
    )
    
    # parameters to be passed in the model
    p = Param(
        # time of simulation
        t_max      = 0.1,
        # directory to save data
        out_dir    = "",
        #CFL
        cfl        = 1.0/16.0,
        # KO diss
        sigma      = 0.0,
        # constraint violation
        damping    = 0.0, # 1.0 (there is constraint damping), or 0.0 (no damping)
        κ1         = 0.0,
        κ2         = 0.0,
        # convention: 1/Mp^2 = 1.0 or = 8*π
        # the 8pi convention does not work well with backreaction; related to CC reg?
        overMp2    = oMp2, #1.0, #8.0*π,
        # cosmological constant; non-zero for backreaction in quantum case
        CC         = CC_theory, #0.0,
        # for Gaussian
        amp        = 0.0, # Minkowski is 0.0
        width      = 2.0,
        rc         = 5.0,
        # infalling_rmax
        infalling_rmax = false,
        # exit the code if an Apparent horizon is found
        AH = false,
        # how often to save data
        save_data   = true,
        data_every  = 8*2^D,
        # how often to save data
        save_data_r0   = false,
        data_r0_every  = 1*2^D,
        # how often to save data for checkpoint
        save_checkpoint  = true,
        checkpoint_every = 1.0, # this is given in hours
        ##########################################################
        # quantum
        hbar = hbar, #1.0, # default = 1.0
        steepness = 1.5,
        r_cut = 20.0,
        # number of quantum modes
        kmax = kmax,
        lmax = lmax,
        # for quantum modes
        dk = dk, #π/15.0, # for the k in the Bessel functions for the quantum ID
        # chose quantum version (regularized or non-regularized)
        PV_reg = true, # false=non-reg., true=reg.
        # masses for PV regularization are in principle [m0, m1, m2, m3, m4, m5]
        # m0 = 0.0
        # m1 = m3, m2 = m4 = sqrt(3)*m1, m5 = 2.0*m1
        # m1 = mPV = 1.0 as default. It can be changed in the example
        # to avoid repetition we have mlist = [m0, m1, m2, m5]
        mlist = [0.0, 1.0*mPV, sqrt(3.0)*mPV, 2.0*mPV],
        # backreaction
        backreaction = false,
        # how often to save bilinears and correlators
        save_quantum      = true,
        quantum_every     = 8*2^D,
        save_quantum_r0   = false,
        quantum_r0_every  = 32*2^D,
        save_bilinears    = true,
        bilinears_every   = 8*2^D,
        save_correlators  = false,
        correlators_every = 16*2^D
    )

    # print messages
    println("kmax = ", p.kmax)
    
    sys = System(g)
    rr = sys.r
    v_classic = zeros(Float64, ( length(rr), 18) )
    # quantum state vector:
    if p.PV_reg==false
        # non-regularized version: one for each k in [1,kmax], l in [0,lmax], 1 + 2 reduction vars 
        v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), 3))
    else
        # regularized version: one for each k in [1,kmax], l in [0,lmax], and m in p.mlist
        v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), length(p.mlist), 3))
    end
    
    println(p.mlist)
    # classic initial data
    v_classic = classical_ID(v_classic, sys, p)
    # quantum initial data
    v_quantum = quantum_ID(v_quantum, sys, p)
    bln = zeros(ComplexF64, (length(rr), 5))
    # calculate the bilinears
    bilinears(0.0, v_classic, v_quantum, p, rr, bln)

    Π  = v_classic[:,2]
    Ψ  = v_classic[:,3]
    A  = v_classic[:,4]
    oA = 1.0./A
    B  = v_classic[:,5]
    oB = 1.0./B
    KB = v_classic[:,9]
    α  = v_classic[:,11]
 
    # classical
    ρ  = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)
    jA = -(oA.^0.5).*oB.*Π.*Ψ
    SA = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)
    SB = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .- Ψ.^2.0)

    ρ_quantum = @. (hbar*c^2/(4.0*π))*((0.5/α^2)*bln[:,2] + (0.5/A)*bln[:,3] + (1/B)*bln[:,5] + 0.5*bln[:,1] );
    jA_quantum = @. -(hbar*c^2/(4.0*π))*(1.0/α)*bln[:,4];
    SA_quantum = @. (hbar*c^2/(4.0*π))*( (0.5/α^2)*bln[:,2] + (0.5/A)*bln[:,3] - (1/B)*bln[:,5] - 0.5*bln[:,1] );
    SB_quantum = @. (hbar*c^2/(4.0*π))*( (0.5/α^2)*bln[:,2] -(0.5/A)*bln[:,3] - 0.5*bln[:,1] );

    ρ_all = (ρ + real(ρ_quantum));
    SA_all =(SA + real(SA_quantum));
    SB_all = (SB + real(SB_quantum));
    jA_all = jA + real.(jA_quantum);

    id = ones(length(rr))
    source_K = @.  (SA_all + 2.0*SB_all - 2.0*CC_theory*id + ρ_all)*oMp2
    source_KB = @. (- SA_all + 2.0*CC_theory*id + ρ_all)*oMp2

    # return the values at r=0 (r[3])
    return source_K[3], source_KB[3] #v_classic, v_quantum, bln, ρ_all, SA_all, SB_all, jA_all, rr

end

kmax_list = collect(20:20:400)
h5write(out_dir*"/kmax_list.h5", "list", kmax_list)

mPV1_source_K = zeros(length(kmax_list))
mPV1_source_KB = zeros(length(kmax_list))

mPV2_source_K = zeros(length(kmax_list))
mPV2_source_KB = zeros(length(kmax_list))

mPV5_source_K = zeros(length(kmax_list))
mPV5_source_KB = zeros(length(kmax_list))

mPV10_source_K = zeros(length(kmax_list))
mPV10_source_KB = zeros(length(kmax_list))

#mPV100_source_K = zeros(length(kmax_list))
#mPV100_source_KB = zeros(length(kmax_list))

@threads for i in 1:length(kmax_list)
    # the syntax (D, lmax, kmax, mPV):
    mPV1_source_K[i], mPV1_source_KB[i] = calculate_bln(3, 10, kmax_list[i], 1);
    mPV2_source_K[i], mPV2_source_KB[i] = calculate_bln(3, 10, kmax_list[i], 2);
    mPV5_source_K[i], mPV5_source_KB[i] = calculate_bln(3, 10, kmax_list[i], 5);
    mPV10_source_K[i], mPV10_source_KB[i] = calculate_bln(3, 10, kmax_list[i], 10);
    #mPV100_source_K[i], mPV100_source_KB[i] = calculate_bln(3, 10, kmax_list[i], 100);
    
end

h5write(out_dir*"/mPV1_source_K.h5", "list", mPV1_source_K)
h5write(out_dir*"/mPV1_source_KB.h5", "list", mPV1_source_KB)

h5write(out_dir*"/mPV2_source_K.h5", "list", mPV2_source_K)
h5write(out_dir*"/mPV2_source_KB.h5", "list", mPV2_source_KB)

h5write(out_dir*"/mPV5_source_K.h5", "list", mPV5_source_K)
h5write(out_dir*"/mPV5_source_KB.h5", "list", mPV5_source_KB)

h5write(out_dir*"/mPV10_source_K.h5", "list", mPV10_source_K)
h5write(out_dir*"/mPV10_source_KB.h5", "list", mPV10_source_KB)

#h5write(out_dir*"/mPV100_source_K.h5", "list", mPV100_source_K)
#h5write(out_dir*"/mPV100_source_KB.h5", "list", mPV100_source_KB)


