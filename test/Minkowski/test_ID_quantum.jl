using Test
using SpheriCo
using Parameters
using SpecialFunctions

# quantum fields initial data
# quantum number k,l, and m is the mass of the scalar field 
function test_uq_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        (k^(l+1)/sqrt(ω))/(gamma(1.5 + l)*2.0^(1+l))
    else
        (k/sqrt(π*ω))*sphericalbesselj(Int(l), k*r)/r^l
    end
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
function test_πq_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        -((sqrt(ω)*k^(l+1))/(gamma(1.5 + l)*2.0^(1+l)))im
    else
        -(k*sqrt(ω/π)*sphericalbesselj(Int(l), k*r)/r^l)im
    end
end

@testset "Test Minkowski ID: quantum" begin

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
    p = SpheriCo.quantum.Param(
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
    dt   = p.cfl* sys.hr
    t = 0.0
    # radial grid
    rr = sys.r

    # without ghost fields
    v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), 3))
    v_quantum = SpheriCo.quantum.quantum_ID(v_quantum, sys, p)
    test_v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), 3))
    for ll in 0:p.lmax
        for kk in 1:p.kmax
            # index for uq, ψq, and πq
            uqi = 1
            ψqi = 2
            πqi = 3
            mi = 0.0
            for i in 3:Nr
                k=kk*p.dk
                # uq; even
                test_v_quantum[i, Int(ll+1), Int(kk), uqi] = test_uq_ID(sys.r[i], k, mi, Float64(ll) )
                # ψq; odd
                test_v_quantum[i, Int(ll+1), Int(kk), ψqi] = test_ψq_ID(sys.r[i], k, mi, Float64(ll) )
                # πq; even
                test_v_quantum[i, Int(ll+1), Int(kk), πqi] = test_πq_ID(sys.r[i], k, mi, Float64(ll) )
            end
            # uq ghosts; even
            test_v_quantum[:, Int(ll+1), Int(kk), uqi] = even_ghosts(test_v_quantum[:, Int(ll+1), Int(kk), uqi])
            # ψq ghosts; odd
            test_v_quantum[:, Int(ll+1), Int(kk), ψqi] = odd_ghosts(test_v_quantum[:, Int(ll+1), Int(kk), ψqi])
            # πq ghosts; even
            test_v_quantum[:, Int(ll+1), Int(kk), πqi] = even_ghosts(test_v_quantum[:, Int(ll+1), Int(kk), πqi])
        end
    end
    # test
    @test maximum(abs.(v_quantum - test_v_quantum)) < 1e-16

    # with ghost fields
    v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), length(p.mlist), 3))
    v_quantum = SpheriCo.quantum.quantum_ID(v_quantum, sys, p)
    test_v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), length(p.mlist), 3))
    for ll in 0:p.lmax
        for kk in 1:p.kmax
            for mm in 1:length(p.mlist)
                # index for uq, ψq, and πq
                uqi = 1
                ψqi = 2
                πqi = 3
                mi = p.mlist[mm]
                for i in 3:Nr
                    k=kk*p.dk
                    # uq; even
                    test_v_quantum[i, Int(ll+1), Int(kk), mm, uqi] = test_uq_ID(sys.r[i], k, mi, Float64(ll) )
                    # ψq; odd
                    test_v_quantum[i, Int(ll+1), Int(kk), mm, ψqi] = test_ψq_ID(sys.r[i], k, mi, Float64(ll) )
                    # πq; even
                    test_v_quantum[i, Int(ll+1), Int(kk), mm, πqi] = test_πq_ID(sys.r[i], k, mi, Float64(ll) )
                end
                # uq ghosts; even
                test_v_quantum[:, Int(ll+1), Int(kk), mm, uqi] = even_ghosts(test_v_quantum[:, Int(ll+1), Int(kk), mm, uqi])
                # ψq ghosts; odd
                test_v_quantum[:, Int(ll+1), Int(kk), mm, ψqi] = odd_ghosts(test_v_quantum[:, Int(ll+1), Int(kk), mm, ψqi])
                # πq ghosts; even
                test_v_quantum[:, Int(ll+1), Int(kk), mm, πqi] = even_ghosts(test_v_quantum[:, Int(ll+1), Int(kk), mm, πqi])
            end
        end
    end
    # test
    @test maximum(abs.(v_quantum - test_v_quantum)) < 1e-16
end