using Test
using SpheriCo
using Parameters
include("../src/classical/parameters.jl")
include("../src/time_integrators.jl")

# the rhs of an advection equation
function F!(t::Float64,
    vf::Array,
    v::Array, 
    sys::System,
    p::Param)

    v = even_ghosts(v)
    # d/dt f = d/dr f
    v_r    = Dr_FD2(v, sys)
    vf[:] = @. -v_r

    nothing
end

# the rhs of an advection equation
# mimics the classical + quantum state vector structure
function F!(t::Float64,
    vc_f::Array,
    vc::Array,
    vq_f::Array,
    vq::Array, 
    sys::System,
    p::Param)

    vc = even_ghosts(vc)
    vq = even_ghosts(vq)
    # d/dt f = d/dr f
    vc_r    = Dr_FD2(vc, sys)
    vc_f[:] = @. -vc_r
    vq_r    = Dr_FD2(vq, sys)
    vq_f[:] = @. -vq_r

    nothing
end

function exact_advect_sol(t::Real, r::Real, amp::Float64, width::Float64, rc::Float64)
    f =  amp*exp(-(((r-t) - rc) / width)^2 ) + amp*exp(-(((r-t) + rc) / width)^2 )
end

@testset "Test time integrators: classic" begin

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
    a = 1.0
    # position of center
    b = 5.0
    #width
    c = 2.0

    # needed to call the RK4 function correctly
    # parameters to be passed in the model
    p = Param(
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

    # initial data
    # "classical"
    vc = zeros(Nr)
    Fc0 = zeros(Nr)
    Fc1 = zeros(Nr)
    Fc2 = zeros(Nr)
    vcf_exact_dt = zeros(Nr) # after time dt
    vcf_exact_3dt = zeros(Nr) # after time dt
    for i in 1:Nr
        vc[i] = SpheriCo.classical.ID_gauss(sys.r[i], a, c, b)
        vcf_exact_dt[i] = exact_advect_sol(dt, sys.r[i], a, c, b)
        vcf_exact_3dt[i] = exact_advect_sol(3*dt, sys.r[i], a, c, b)
    end

    # save the rhs at t=0
    F!(t, Fc0, vc, sys, p)
    # one timestep forward by dt; changes vf
    RK4(t, dt, vc, sys, p)
    t = t+dt

    # test the RK4 sol at t=dt
    @test maximum(abs.(vc - vcf_exact_dt)[3:end]) < 1e-9

    # continoue to test the AB3!

    # save the rhs at t=dt
    F!(t, Fc1, vc, sys, p)
    # one timestep forward by dt; changes vf
    RK4(t, dt, vc, sys, p)
    t = t+dt
    # save the rhs at t=2*dt
    F!(t, Fc2, vc, sys, p)
    
    # get the AB3! sol at t=3*dt
    AB3!(vc, Fc2, Fc1, Fc0, dt)

    @test maximum(abs.(vc - vcf_exact_3dt)[3:end]) < 1e-8

end

@testset "Test time integrators: classical + quantum" begin

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
    a = 1.0
    # position of center
    b = 5.0
    #width
    c = 2.0

    # needed to call the RK4 function correctly
    # parameters to be passed in the model
    p = Param(
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

    # initial data
    # "classical"
    vc = zeros(Nr)
    Fc0 = zeros(Nr)
    Fc1 = zeros(Nr)
    Fc2 = zeros(Nr)
    vcf_exact_dt = zeros(Nr) # after time dt
    vcf_exact_3dt = zeros(Nr) # after time dt
    # "quantum"
    vq = zeros(Nr)
    Fq0 = zeros(Nr)
    Fq1 = zeros(Nr)
    Fq2 = zeros(Nr)
    vqf_exact_dt = zeros(Nr) # after time dt
    vqf_exact_3dt = zeros(Nr) # after time dt
    for i in 1:Nr
        vc[i] = SpheriCo.classical.ID_gauss(sys.r[i], a, c, b)
        vcf_exact_dt[i] = exact_advect_sol(dt, sys.r[i], a, c, b)
        vcf_exact_3dt[i] = exact_advect_sol(3*dt, sys.r[i], a, c, b)
        vq[i] = SpheriCo.classical.ID_gauss(sys.r[i], a, c, b)
        vqf_exact_dt[i] = exact_advect_sol(dt, sys.r[i], a, c, b)
        vqf_exact_3dt[i] = exact_advect_sol(3*dt, sys.r[i], a, c, b)
    end

    # save the rhs at t=0
    F!(t, Fc0, vc, Fq0, vq, sys, p)
    # one timestep forward by dt; changes vf
    RK4(t, dt, vc, vq, sys, p)
    t = t+dt

    # test the RK4 sol at t=dt
    @test maximum(abs.(vc - vcf_exact_dt)[3:end]) < 1e-9
    @test maximum(abs.(vq - vqf_exact_dt)[3:end]) < 1e-9

    # continoue to test the AB3!

    # save the rhs at t=dt
    F!(t, Fc1, vc, Fq1, vq, sys, p)
    # one timestep forward by dt; changes vf
    RK4(t, dt, vc, vq, sys, p)
    t = t+dt
    # save the rhs at t=2*dt
    F!(t, Fc2, vc, Fq2, vq, sys, p)
    
    # get the AB3! sol at t=3*dt
    AB3!(vc, Fc2, Fc1, Fc0, vq, Fq2, Fq1, Fq0, dt)

    @test maximum(abs.(vc - vcf_exact_3dt)[3:end]) < 1e-8
    @test maximum(abs.(vq - vqf_exact_3dt)[3:end]) < 1e-8

end
