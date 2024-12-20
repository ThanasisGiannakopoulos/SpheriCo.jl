using Test
using SpheriCo
using Parameters
using Interpolations

@testset "Test classical rhs: Minkowski" begin

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
    p = SpheriCo.classical.Param(
    # time of simulation
    t_max      = 1,
    # directory to save data
    out_dir    = "./",
    #CFL
    cfl        = 1.0/8,
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
    # random data amplitude (for robust stability test)
    rand = false,
    A_rand = 1.0,
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
    checkpoint_every = 1.0 # this is given in hours
    )

    # timestep
    dt   = p.cfl* sys.hr
    t = 0.0
    # radial grid
    rr = sys.r
    # v_classic = [Φ, Π, Ψ, A, B, DB, Utld, K, KB, λ, α, Dα, Θ, Zr, f, g, U, V]^T
    v_classic = zeros(Float64, ( length(rr), 18) )
    v_classic = SpheriCo.classical.classical_ID(v_classic, sys, p)
    # initiate the rhs as random numbers
    rhs = rand(Float64, ( length(rr), 18) )
    SpheriCo.classical.F!(t, rhs, v_classic, sys, p)
    # test rhs
    @test maximum(abs.(rhs[:,1:16] - zeros(Float64, (length(rr), 16)))[3:end]) < 1e-16
    @test maximum(abs.(rhs[:,17:18] - ones(Float64, (length(rr), 2)))[3:end]) < 1e-16
end