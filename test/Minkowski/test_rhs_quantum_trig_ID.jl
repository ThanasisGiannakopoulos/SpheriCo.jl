using Test
using SpheriCo
using Parameters
using Interpolations

# without ghosts
@testset "Test quantum rhs: Minkowski (no ghosts, no backreaction, cos/sin quantum ID)" begin

    # change D for number of points
    D = 7
    Nr = (128)*2^D + 1 + 2 # +1 for last point, + 2 for ghosts left of r=0
    # positions of outer boundary
    rmax  = 10
    # radial grid
    g = Grid(
    # discretization parameters
    Nr         = Nr,
    r_max      = rmax
    )

    # pass the parameters of the system
    sys = System(g)

    # gaussian amplitude
    a = 0.0
    # position of center
    b = 0.0
    #width
    c = 2.0

    # parameters to be passed in the model
    pc = SpheriCo.classical.Param(
    # time of simulation
    t_max      = 1,
    # directory to save data
    out_dir    = "./",
    #CFL
    cfl        = 1.0/16,
    # KO diss
    sigma      = 0.0,
    # constraint violation
    damping    = 0.0, # 1.0 (there is constraint damping), or 0.0 (no damping)
    κ1         = 0.0,
    κ2         = 0.0,
    # convention: 1/Mp^2 = 1.0 or = 8*π
    overMp2    = 1.0, # 8.0*π or 1.0,
    # cosmological constant; non-zero for backreaction in quantum case
    CC         = 0.0,
    # for Gaussian
    amp        = a,
    rc         = b,
    width      = c,
    # infalling_rmax
    infalling_rmax = false,
    # to exit the code if an Apparent horizon is found
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
    )

    # parameters to be passed in the model
    pq = SpheriCo.quantum.Param(
    # time of simulation
    t_max      = 1,
    # directory to save data
    out_dir    = "./",
    #CFL
    cfl        = 1.0/16,
    # KO diss
    sigma      = 0.0,
    # constraint violation
    damping    = 0.0, # 1.0 (there is constraint damping), or 0.0 (no damping)
    κ1         = 0.0,
    κ2         = 0.0,
    # convention: 1/Mp^2 = 1.0 or = 8*π
    overMp2    = 1.0, # 8.0*π or 1.0,
    # cosmological constant; non-zero for backreaction in quantum case
    CC         = 0.0,
    # for Gaussian
    amp        = a,
    rc         = b,
    width      = c,
    # infalling_rmax
    infalling_rmax = false,
    # to exit the code if an Apparent horizon is found
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
    hbar = 1.0,
    # for filter based on tanh, to cut backreaction in causally disconnected r-part
    # needed for stability
    steepness = 1.0,
    # also for filter; the r where the transition from 1->0 happens
    r_cut = 20.0,
    # number of quantum modes
    kmax = 5.0,
    lmax = 5.0,
    # for quantum modes
    dk = π/10.0, # for the k in the Bessel functions for the quantum ID
    # chose quantum version (regularized or non-regularized)
    PV_reg = false, # false -> non-reg., true -> reg.
    # masses for PV regularization are in principle [m0, m1, m2, m3, m4, m5]
    mlist = [0.0, 1.0, sqrt(3.0), 2.0],
    # backreaction
    backreaction = false,
    # how often to save bilinears and correlators
    save_quantum = true,
    quantum_every = 1,
    save_quantum_r0 = false,
    quantum_r0_every = 1,
    save_bilinears = false,
    bilinears_every = 1,
    save_correlators = false,
    correlators_every = 1,
    )

    # timestep
    dt   = pq.cfl* sys.hr
    t = 0.0
    # radial grid
    rr = sys.r
    # iniial data
    # v_classic = [Φ, Π, Ψ, A, B, DB, Utld, K, KB, λ, α, Dα, Θ, Zr, f, g, U, V]^T
    v_classic = zeros(Float64, ( length(rr), 18) )
    v_classic = SpheriCo.classical.classical_ID(v_classic, sys, pc)
    v_quantum = zeros(ComplexF64, (length(rr), Int(pq.lmax+1), Int(pq.kmax), 3))
    for kk in 1:Int(pq.kmax)
        for l in 0:Int(pq.lmax)
            # uq ID
            v_quantum[:,l+1,kk,1] = cos.(rr) + 0.0im*rr
            # ψq ID
            v_quantum[:,l+1,kk,2] = -sin.(rr) + 0.0im*rr
            # πq ID
            v_quantum[:,l+1,kk,3] = 0.0*rr - cos.(rr)*1.0im
        end
    end
    # initiate the rhs as random numbers
    rhs_classic = rand(Float64, (length(rr), 18) )
    rhs_quantum = rand(ComplexF64, (length(rr), Int(pq.lmax+1), Int(pq.kmax), 3))
    SpheriCo.quantum.F!(t, rhs_classic, v_classic, rhs_quantum, v_quantum, sys, pq)
    # test rhs_classic
    @test maximum(abs.(rhs_classic[:,1:16] - zeros(Float64, (length(rr), 16)))[3:end]) < 1e-16
    @test maximum(abs.(rhs_classic[:,17:18] - ones(Float64, (length(rr), 2)))[3:end]) < 1e-16
    # test rhs_quantum
    # mass m = 0
    m = 0.0
    for k in 1:Int(pq.kmax)
        kk = pq.dk*k
        for l in 0:Int(pq.lmax)
            ll = float(l)
            # uq rhs = πq
            @test maximum(abs.(rhs_quantum[:,l+1,k,1] - (0.0*rr - cos.(rr)im))[3:end]) < 1e-16
            # ψq rhs = πq'
            @test maximum(abs.(rhs_quantum[:,l+1,k,2] - (0.0*rr + sin.(rr)*1.0im))[3:end-1]) < 1e-7
            # πq rhs = 2(l+1)*ψq/r + ψq' - m^2*uq ; m=0 here
            @test maximum(abs.(rhs_quantum[:,l+1,k,3][4:end-1] - 2*(ll+1)*(-sin.(rr) + 0.0im*rr)[4:end-1]./rr[4:end-1] - (-cos.(rr) + 0.0im*rr)[4:end-1])) < 1e-5
            # at r=0
            @test maximum(abs.(rhs_quantum[:,l+1,k,3][3] - 2*(ll+1+0.5)*(-cos(0) + 0.0im))) < 1e-5
        end
    end
end

# with ghosts
@testset "Test quantum rhs: Minkowski (no backreaction, cos/sin quantum ID)" begin

    # change D for number of points
    D = 7
    Nr = (128)*2^D + 1 + 2 # +1 for last point, + 2 for ghosts left of r=0
    # positions of outer boundary
    rmax  = 10
    # radial grid
    g = Grid(
    # discretization parameters
    Nr         = Nr,
    r_max      = rmax
    )

    # pass the parameters of the system
    sys = System(g)

    # gaussian amplitude
    a = 0.0
    # position of center
    b = 0.0
    #width
    c = 2.0

    # parameters to be passed in the model
    pc = SpheriCo.classical.Param(
    # time of simulation
    t_max      = 1,
    # directory to save data
    out_dir    = "./",
    #CFL
    cfl        = 1.0/16,
    # KO diss
    sigma      = 0.0,
    # constraint violation
    damping    = 0.0, # 1.0 (there is constraint damping), or 0.0 (no damping)
    κ1         = 0.0,
    κ2         = 0.0,
    # convention: 1/Mp^2 = 1.0 or = 8*π
    overMp2    = 1.0, # 8.0*π or 1.0,
    # cosmological constant; non-zero for backreaction in quantum case
    CC         = 0.0,
    # for Gaussian
    amp        = a,
    rc         = b,
    width      = c,
    # infalling_rmax
    infalling_rmax = false,
    # to exit the code if an Apparent horizon is found
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
    )

    # parameters to be passed in the model
    pq = SpheriCo.quantum.Param(
    # time of simulation
    t_max      = 1,
    # directory to save data
    out_dir    = "./",
    #CFL
    cfl        = 1.0/16,
    # KO diss
    sigma      = 0.0,
    # constraint violation
    damping    = 0.0, # 1.0 (there is constraint damping), or 0.0 (no damping)
    κ1         = 0.0,
    κ2         = 0.0,
    # convention: 1/Mp^2 = 1.0 or = 8*π
    overMp2    = 1.0, # 8.0*π or 1.0,
    # cosmological constant; non-zero for backreaction in quantum case
    CC         = 0.0,
    # for Gaussian
    amp        = a,
    rc         = b,
    width      = c,
    # infalling_rmax
    infalling_rmax = false,
    # to exit the code if an Apparent horizon is found
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
    hbar = 1.0,
    # for filter based on tanh, to cut backreaction in causally disconnected r-part
    # needed for stability
    steepness = 1.0,
    # also for filter; the r where the transition from 1->0 happens
    r_cut = 20.0,
    # number of quantum modes
    kmax = 5.0,
    lmax = 5.0,
    # for quantum modes
    dk = π/10.0, # for the k in the Bessel functions for the quantum ID
    # chose quantum version (regularized or non-regularized)
    PV_reg = true, # false -> non-reg., true -> reg.
    # masses for PV regularization are in principle [m0, m1, m2, m3, m4, m5]
    mlist = [0.0, 1.0, sqrt(3.0), 2.0],
    # backreaction
    backreaction = false,
    # how often to save bilinears and correlators
    save_quantum = true,
    quantum_every = 1,
    save_quantum_r0 = false,
    quantum_r0_every = 1,
    save_bilinears = false,
    bilinears_every = 1,
    save_correlators = false,
    correlators_every = 1,
    )

    # timestep
    dt   = pq.cfl* sys.hr
    t = 0.0
    # radial grid
    rr = sys.r
    # iniial data
    # v_classic = [Φ, Π, Ψ, A, B, DB, Utld, K, KB, λ, α, Dα, Θ, Zr, f, g, U, V]^T
    v_classic = zeros(Float64, ( length(rr), 18) )
    v_classic = SpheriCo.classical.classical_ID(v_classic, sys, pc)
    v_quantum = zeros(ComplexF64, (length(rr), Int(pq.lmax+1), Int(pq.kmax), 4, 3))
    for kk in 1:Int(pq.kmax)
        for l in 0:Int(pq.lmax)
            for mm in 1:4
                # uq ID
                v_quantum[:,l+1,kk,mm,1] = cos.(rr) + 0.0im*rr
                # ψq ID
                v_quantum[:,l+1,kk,mm,2] = -sin.(rr) + 0.0im*rr
                # πq ID
                v_quantum[:,l+1,kk,mm,3] = 0.0*rr - cos.(rr)*1.0im
            end
        end
    end
    # initiate the rhs as random numbers
    rhs_classic = rand(Float64, (length(rr), 18) )
    rhs_quantum = rand(ComplexF64, (length(rr), Int(pq.lmax+1), Int(pq.kmax), 4, 3))
    SpheriCo.quantum.F!(t, rhs_classic, v_classic, rhs_quantum, v_quantum, sys, pq)
    # test rhs_classic
    @test maximum(abs.(rhs_classic[:,1:16] - zeros(Float64, (length(rr), 16)))[3:end]) < 1e-16
    @test maximum(abs.(rhs_classic[:,17:18] - ones(Float64, (length(rr), 2)))[3:end]) < 1e-16
    # test rhs_quantum
    # skip point at rmax due to outgoing boundary condition
    for k in 1:Int(pq.kmax)
        kk = pq.dk*k
        for l in 0:Int(pq.lmax)
            for mm in 1:4
                m = pq.mlist[mm]
                ll = float(l)
                # uq rhs = πq
                @test maximum(abs.(rhs_quantum[:,l+1,k,mm,1] - (0.0*rr - cos.(rr)im))[3:end-1]) < 1e-16
                # ψq rhs = πq'
                @test maximum(abs.(rhs_quantum[:,l+1,k,mm,2] - (0.0*rr + sin.(rr)*1.0im))[3:end-1]) < 1e-7
                # πq rhs = 2(l+1)*ψq/r + ψq' - m^2*uq ; m=0 here
                @test maximum(abs.(rhs_quantum[:,l+1,k,mm,3][4:end-1] - 2*(ll+1)*(-sin.(rr) + 0.0im*rr)[4:end-1]./rr[4:end-1] - (-cos.(rr) + 0.0im*rr)[4:end-1] + m^2*(cos.(rr) + 0.0im*rr)[4:end-1])) < 1e-5
                # at r=0
                @test maximum(abs.(rhs_quantum[:,l+1,k,mm,3][3] - 2*(ll+1+0.5)*(-cos(0) + 0.0im) + m^2*(cos(0) + 0.0im))) < 1e-5
            end
        end
    end
end