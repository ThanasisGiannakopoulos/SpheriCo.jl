using Test
using SpheriCo
using Parameters
using Interpolations
using SpecialFunctions

function test_uq_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        (k^(l+1)/sqrt(ω))/(gamma(1.5 + l)*2.0^(1+l))
    else
        (k/sqrt(π*ω))*sphericalbesselj(Int(l), k*r)/r^l
    end
end
function test_uq_ID(r::Array, k::Float64, m::Float64, l::Float64)
    f = zeros(ComplexF64, length(r))
    for i in 3:length(r)
        f[i] = test_uq_ID(r[i], k, m, l)
    end
    return f
end

function test_ψq_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        0.0 + 0.0im
    else
        (k/sqrt(π*ω))*0.5*(
            k*r*sphericalbesselj(Int(l)-1,k*r) -
            (2*l+1)*sphericalbesselj(Int(l),k*r) -
            k*r*sphericalbesselj(Int(l)+1,k*r) 
        )/(r^(l+1))
    end
end
function test_ψq_ID(r::Array, k::Float64, m::Float64, l::Float64)
    f = zeros(ComplexF64, length(r))
    for i in 3:length(r)
        f[i] = test_ψq_ID(r[i], k, m, l)
    end
    return f
end

function test_πq_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        -((sqrt(ω)*k^(l+1))/(gamma(1.5 + l)*2.0^(1+l)))im
    else
        -(k*sqrt(ω/π)*sphericalbesselj(Int(l), k*r)/r^l)im
    end
end
function test_πq_ID(r::Array, k::Float64, m::Float64, l::Float64)
    f = zeros(ComplexF64, length(r))
    for i in 3:length(r)
        f[i] = test_πq_ID(r[i], k, m, l)
    end
    return f
end

function test_πq_r_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        0.0 + 0.0im
    else
        -((k*sqrt(ω)/sqrt(π))*0.5*(
            k*r*sphericalbesselj(Int(l)-1,k*r) -
            (2*l+1)*sphericalbesselj(Int(l),k*r) -
            k*r*sphericalbesselj(Int(l)+1,k*r) 
        )/(r^(l+1)))im
    end
end
function test_πq_r_ID(r::Array, k::Float64, m::Float64, l::Float64)
    f = zeros(ComplexF64, length(r))
    for i in 3:length(r)
        f[i] = test_πq_r_ID(r[i], k, m, l)
    end
    return f
end

function test_ψq_r_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        0.0 + 0.0im # this is wrong
    else
        (k/sqrt(π*ω))*(
            -k*l*r*sphericalbesselj(Int(l)-1,k*r) +
            (l+2*l^2-k^2*r^2)*sphericalbesselj(Int(l),k*r) +
            k*(2+l)*r*sphericalbesselj(Int(l)+1,k*r)
        )/(r^(l+2))
    end
end
function test_ψq_r_ID(r::Array, k::Float64, m::Float64, l::Float64)
    f = zeros(ComplexF64, length(r))
    for i in 3:length(r)
        f[i] = test_ψq_r_ID(r[i], k, m, l)
    end
    return f
end

# without ghosts
@testset "Test quantum rhs: Minkowski (no ghosts, no backreaction)" begin

    # change D for number of points
    D = 6
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
    v_quantum = SpheriCo.quantum.quantum_ID(v_quantum, sys, pq)

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
            @test maximum(abs.(rhs_quantum[:,l+1,k,1] - test_πq_ID(rr,kk,m,ll))[3:end]) < 1e-16
            # ψq rhs = πq'
            @test maximum(abs.(rhs_quantum[:,l+1,k,2] - test_πq_r_ID(rr,kk,m,ll))[3:end]) < 2*1e-5
            # πq rhs = 2(l+1)*ψq/r - ψq' - m^2*uq ; m=0 here
            # at r=0 the test is wrong
            @test maximum(abs.(rhs_quantum[:,l+1,k,3][4:end] - 2*(ll+1)*test_ψq_ID(rr,kk,m,ll)[4:end]./rr[4:end] - test_ψq_r_ID(rr,kk,m,ll)[4:end])) < 2*1e-4
        end
    end
end