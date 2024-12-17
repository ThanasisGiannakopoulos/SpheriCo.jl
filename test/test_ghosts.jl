using Test
using SpheriCo

@testset "ghost tests" begin

    # change D for number of points
    D = 4
    Nr = (128)*2^D + 1 + 2 # +1 for last point, + 2 for ghosts left of r=0
    # positions of outer boundary
    rmax  = 10
    # radial grid
    g = Grid(
    # discretization parameters
    Nr         = Nr,
    r_max      = rmax
    )
    
    # generate random data
    f = rand(Nr)
    g = rand(Nr)

    # populate ghosts on the left of r=0
    f = even_ghosts(f)
    g = odd_ghosts(g)

    @test abs(f[4]-f[2]) < 1e-16
    @test abs(f[5]-f[1]) < 1e-16
    @test abs(g[4]+g[2]) < 1e-16
    @test abs(g[5]+g[1]) < 1e-16
end
