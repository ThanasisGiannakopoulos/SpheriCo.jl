
using SpheriCo
using SpheriCo.classical

"""
The variables below are global and fed later into the paramterersof
the model.  They are used to build the name of the folder where the
simulation is saved.  A different choice is to manually give the name
of that folder, and tune the paramteres below, directly.
"""

# amplitude
"""
relation for different conventions:
if overMp2 = 8.0*π  and amp = A
then the same is 
overMp2 = 1 and amp = A*sqrt(8.0*π)
"""
a = 0.25
# position of center
b = 0.0
#width
c = 4.0

# change D for number of points
D = 3
Nr = (128)*2^D + 1 + 2 # +1 for last point, + 2 for ghosts left of r=0
# for artificial dissipation
sigma = 0.0
# for timestep dt=dr/cfl_denom
cfl_denom = 16.0
# for constraint damping; 0.0 is off and 1.0 is on
damp = 0.0
# control the amount of constraint damping; need to be >1
kk1 = 0.0
kk2 = 0.0
# positions of outer boundary
rmax  = 15
# time of the evolution
tmax  = 7.0
# simulate with an infalling r_max (at the speed of light); true or false
infalling_rmax = false
# used to bread
AH = true
# for random noise tests
rand = false
A_rand_base = 10^(-3)
A_rand = A_rand_base/4^D
# convention for 1/M_p^2
overMp2 = 8.0*π

# name of the folder where the data are saved
if rand==true && damp==0.0
    root_dir = "./classical_runs/a$(a)_b$(b)_c$(c)_rmax$(rmax)_tmax$(tmax)_cfl$(1.0/cfl_denom)_sigma$(sigma)_infalling_rmax_$(infalling_rmax)_rand_$(rand)_$(A_rand_base)_overMp2_$(overMp2)_damp$(Int(damp))/"
end
if rand==true && damp==1.0
    root_dir = "./classical_runs/a$(a)_b$(b)_c$(c)_rmax$(rmax)_tmax$(tmax)_cfl$(1.0/cfl_denom)_sigma$(sigma)_infalling_rmax_$(infalling_rmax)_rand_$(rand)_$(A_rand_base)_overMp2_$(overMp2)_damp$(Int(damp))_k1_$(kk1)_k2_$(kk2)/"
end
if rand==false && damp==0.0
    root_dir = "./classical_runs/a$(a)_b$(b)_c$(c)_rmax$(rmax)_tmax$(tmax)_cfl$(1.0/cfl_denom)_sigma$(sigma)_infalling_rmax_$(infalling_rmax)_rand_$(rand)_overMp2_$(overMp2)_damp$(Int(damp))/"
end
if rand==false && damp==1.0
    root_dir = "./classical_runs/a$(a)_b$(b)_c$(c)_rmax$(rmax)_tmax$(tmax)_cfl$(1.0/cfl_denom)_sigma$(sigma)_infalling_rmax_$(infalling_rmax)_rand_$(rand)_overMp2_$(overMp2)_damp$(Int(damp))_k1_$(kk1)_k2_$(kk2)/"
end

# create the folders where data are saved
out_dir = joinpath(root_dir, "data_$(Nr)/")
mkpath(out_dir)
cp("./evol_classical.jl", out_dir*"/evol_classical.jl", force=true)

# radial grid
g = Grid(
    # discretization parameters
    Nr         = Nr,
    r_max      = rmax
)

# parameters to be passed in the model
p = Param(
    # time of simulation
    t_max      = tmax,
    # directory to save data
    out_dir    = out_dir,
    #CFL
    cfl        = 1.0/cfl_denom,
    # KO diss
    sigma      = sigma,
    # constraint violation
    damping    = damp, # 1.0 (there is constraint damping), or 0.0 (no damping)
    κ1         = kk1,
    κ2         = kk2,
    # convention: 1/Mp^2 = 1.0 or = 8*π
    overMp2    = overMp2, # 8.0*π or 1.0,
    # cosmological constant; non-zero for backreaction in quantum case
    CC         = 0.0,
    # for Gaussian
    amp        = a,
    rc         = b,
    width      = c,
    # random data amplitude (for robust stability test)
    rand = rand,
    A_rand = A_rand,
    # infalling_rmax
    infalling_rmax = infalling_rmax,
    # to exit the code if an Apparent horizon is found
    AH = AH,
    # how often to save data
    save_data   = true,
    data_every  = 2*2^D,
    # how often to save data
    save_data_r0   = false,
    data_r0_every  = 1*2^D,
    # how often to save data for checkpoint
    save_checkpoint  = true,
    checkpoint_every = 1.0 # this is given in hours
)

run_classical(g, p)
