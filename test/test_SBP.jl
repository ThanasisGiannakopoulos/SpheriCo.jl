using Test
using SpheriCo

"""
test the approximation (SBP2l-FD2)@f == l*f/r; with f an odd function
TODO: it is not checking the last point at rmax (there is an error there maybe)
practically we ignore the outer part of the domain in physical interpretation, but it could be improved
"""

# define the 2nd radial derivative of the gaussian IDfunction ID_gauss_r(r::Real, amp::Float64, width::Float64, rc::Float64)
function ID_gauss_rr(r::Real, amp::Float64, width::Float64, rc::Float64)
    f = (2*amp*exp(-((r + rc)^2/width^2))*(2*(1 + exp((4*r*rc)/width^2))*r^2 -
     4*(-1 + exp((4*r*rc)/width^2))*r*rc +
     (1 + exp((4*r*rc)/width^2))*(2*rc^2 - width^2)))/width^4
end

@testset "SBP tests (real valued)" begin
    # gaussian amplitude
    a = 1.0
    # position of center
    b = 0.0
    #width
    c = 2.0
    
    # change D for number of points
    D = 9
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
    # initiate the test function to d/dr of Gaussian
    f = zeros(Nr)
    for i in 1:Nr
        f[i] = SpheriCo.classical.ID_gauss_r(sys.r[i], a, c, b)
    end
    f = odd_ghosts(f)

    # l = 1
    l = 1
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-6

    # l = 10
    l = 10
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-6

    # l = 40
    l = 40
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-5

    # l = 80
    l = 80
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-4

    # l = 100
    l = 100
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-4

end

@testset "SBP tests (complex valued)" begin
    # gaussian amplitude
    a = 1.0
    # position of center
    b = 0.0
    #width
    c = 2.0
    
    # change D for number of points
    D = 9
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
    # initiate the test function to d/dr of Gaussian
    f = zeros(ComplexF64, Nr)
    for i in 1:Nr
        f[i] = SpheriCo.classical.ID_gauss_r(sys.r[i], a, c, b) + SpheriCo.classical.ID_gauss_r(sys.r[i], a, c, b)*im
    end
    f = odd_ghosts(f)

    # l = 1
    l = 1
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*(ID_gauss_rr(0,a, c, b) + ID_gauss_rr(0,a, c, b)*im)
    @test maximum(abs.(real(diff[3:end-1]))) < 1e-6
    @test maximum(abs.(imag(diff[3:end-1]))) < 1e-6

    # l = 10
    l = 10
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*(ID_gauss_rr(0,a, c, b) + ID_gauss_rr(0,a, c, b)*im)
    @test maximum(abs.(real(diff[3:end-1]))) < 1e-6
    @test maximum(abs.(imag(diff[3:end-1]))) < 1e-6

    # l = 40
    l = 40
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*(ID_gauss_rr(0,a, c, b) + ID_gauss_rr(0,a, c, b)*im)
    @test maximum(abs.(real(diff[3:end-1]))) < 1e-5
    @test maximum(abs.(imag(diff[3:end-1]))) < 1e-5

    # l = 80
    l = 80
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*(ID_gauss_rr(0,a, c, b) + ID_gauss_rr(0,a, c, b)*im)
    @test maximum(abs.(real(diff[3:end-1]))) < 1e-4
    @test maximum(abs.(imag(diff[3:end-1]))) < 1e-4

    # l = 100
    l = 100
    diff = Dr_SBP2(f, sys, l) - Dr_FD2(f, sys) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*(ID_gauss_rr(0,a, c, b) + ID_gauss_rr(0,a, c, b)*im)
    @test maximum(abs.(real(diff[3:end-1]))) < 1e-4
    @test maximum(abs.(imag(diff[3:end-1]))) < 1e-4

end

@testset "SBP tests (postprocessing version)" begin
    # gaussian amplitude
    a = 1.0
    # position of center
    b = 0.0
    #width
    c = 2.0
    
    # change D for number of points
    D = 9
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
    hr = sys.r[4] - sys.r[3]
    # initiate the test function to d/dr of Gaussian
    f = zeros(Nr)
    for i in 1:Nr
        f[i] = SpheriCo.classical.ID_gauss_r(sys.r[i], a, c, b)
    end
    f = odd_ghosts(f)

    # l = 1
    l = 1
    diff = Dr_SBP2(f, hr, l) - Dr_FD2(f, hr) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, hr, l)[3] - Dr_FD2(f, hr)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-6

    # l = 10
    l = 10
    diff = Dr_SBP2(f, hr, l) - Dr_FD2(f, hr) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, hr, l)[3] - Dr_FD2(f, hr)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-6

    # l = 40
    l = 40
    diff = Dr_SBP2(f, hr, l) - Dr_FD2(f, hr) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, sys, l)[3] - Dr_FD2(f, sys)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-5

    # l = 80
    l = 80
    diff = Dr_SBP2(f, hr, l) - Dr_FD2(f, hr) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, hr, l)[3] - Dr_FD2(f, hr)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-4

    # l = 100
    l = 100
    diff = Dr_SBP2(f, hr, l) - Dr_FD2(f, hr) - l*f./sys.r
    # l'Hospital at r=0
    diff[3] = Dr_SBP2(f, hr, l)[3] - Dr_FD2(f, hr)[3] - l*ID_gauss_rr(0,a, c, b)
    @test maximum(abs.(diff[3:end-1])) < 1e-4

end