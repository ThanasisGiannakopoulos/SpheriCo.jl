using Test
using SpheriCo

# define the 2nd radial derivative of the gaussian ID 
function ID_gauss_rr(r::Real, amp::Float64, width::Float64, rc::Float64)
    f = (2*amp*exp(-((r + rc)^2/width^2))*(2*(1 + exp((4*r*rc)/width^2))*r^2 -
     4*(-1 + exp((4*r*rc)/width^2))*r*rc +
     (1 + exp((4*r*rc)/width^2))*(2*rc^2 - width^2)))/width^4
end

@testset "FD2 test" begin
    # gaussian amplitude
    a = 1.0
    # position of center
    b = 10.0
    #width
    c = 2.0
    
    # change D for number of points
    D = 6
    Nr = (128)*2^D + 1 + 2 # +1 for last point, + 2 for ghosts left of r=0
    # positions of outer boundary
    rmax  = 20
    # radial grid
    g = Grid(
    # discretization parameters
    Nr         = Nr,
    r_max      = rmax
    )
    
    # pass the parameters of the system
    sys = System(g)

    f = zeros(Nr)
    fr_exact = zeros(Nr)
    frr_exact = zeros(Nr)
    for i in 1:Nr
        f[i] = SpheriCo.classical.ID_gauss(sys.r[i], a, c, b)
        fr_exact[i] = SpheriCo.classical.ID_gauss_r(sys.r[i], a, c, b)    
        frr_exact[i] = ID_gauss_rr(sys.r[i], a, c, b)    
    end

    fr1 = Dr_FD2(f, sys) # in the code
    frr = Drr_FD2(f, sys) # in the code
    hr = sys.r[4] - sys.r[3]
    fr2 = Dr_FD2(f, hr) # postprocessing version

    #@test fr ≈ fr_exact
    @test maximum(abs.(fr1 - fr_exact)) < 1e-6  # Tolerance can be adjusted
    @test maximum(abs.(fr2 - fr_exact)) < 1e-6  # Tolerance can be adjusted
    @test maximum(abs.(frr - frr_exact)) < 1e-6  # Tolerance can be adjusted

end
