@with_kw struct Param
    #same as classical
    # max evolution time
    t_max      :: Float64 = 10.0
    # directory to save data
    out_dir    :: String
    #CFL
    cfl        :: Float64 = 0.125
    # KO diss
    sigma      :: Float64 = 0.02
    # constraint violation
    damping    :: Float64 = 0.0 # 1.0 (there is constraint damping), or 0.0 (no damping)
    κ1         :: Float64 = 0.02
    κ2         :: Float64 = 0.0
    # convention: 1/Mp^2 = 1.0 or = 8*π
    overMp2    :: Float64 = 8*π
    # cosmological constant; non-zero for quantum with backreaction
    CC         :: Float64 = 0.0
    # classical Gaussian profile
    amp        :: Float64
    width      :: Float64
    rc         :: Float64
    # infalling_rmax
    infalling_rmax :: Bool = false
    # to exit the code if an Apparent horizon is found
    AH :: Bool = false
    # how often to save data
    save_data   :: Bool = true
    data_every  :: Int
    # how often to save data at r0
    save_data_r0   :: Bool = false
    data_r0_every  :: Int
    # how often to save checkpoint
    save_checkpoint  :: Bool = true
    checkpoint_every :: Float64 = 1.0 # this is given in hours
    ##########################################################
    # quantum
    hbar :: Float64 = 1.0
    # for filter based on tanh, to cut backreaction in causally disconnecter r-part
    # needed for stability
    steepness ::Float64 = 1.5
    # also for filter; the r where the transition from 1->0 happens
    r_cut ::Float64
    # number of quantum modes
    kmax       :: Float64
    lmax       :: Float64
    # for quantum modes
    dk         :: Float64 = π/15.0 # for the k in the Bessel functions for the quantum ID
    # chose quantum version (regularized or non-regularized)
    PV_reg :: Bool = false # false=non-reg., true=reg.
    # masses for PV regularization are in principle [m0, m1, m2, m3, m4, m5]
    # m0 = 0.0
    # m1 = m3, m2 = m4 = sqrt(3)*m1, m5 = 2.0*m1
    # m1 = mPV = 1.0 as default. It can be changed in the example
    # to avoid repetition we have mlist = [m0, m1, m2, m5]
    mlist :: Vector{Float64} = [0.0, 1.0, sqrt(3.0), 2.0]
    # backreaction
    backreaction::Bool = false
    # how often to save bilinears and correlators
    save_quantum      :: Bool = true
    quantum_every     :: Int
    save_quantum_r0   :: Bool = false
    quantum_r0_every  :: Int
    save_bilinears    :: Bool = true
    bilinears_every   :: Int
    save_correlators  :: Bool = true
    correlators_every :: Int
    
end
